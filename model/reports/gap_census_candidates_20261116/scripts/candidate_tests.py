"""LP test of every single-reaction candidate for every root metabolite of candidate 7ccbb568.

Usage:
    python candidate_tests.py test   MODEL.xml ART_DIR REAC_PROP.tsv GAPFILL_DIR
    python candidate_tests.py tables ART_DIR REPORT_DIR DECISIONS.csv

test    checks the input hashes (../input_hashes.csv; needs IYALI26_RESEARCH_ROOT), builds the candidate list
        from the census root table, checks it against the plan-time count, runs the LP tests with GLPK and
        Gurobi and writes ART_DIR/out/candidate_*.{csv,json}
tables  merges the test results with the census fix tests into the report tables; stops when the test stage
        recorded a consistency problem or a census fix test failed to solve

A candidate for a root is a MetaNetX non-transport reaction that contains the root's compound, whose every
compound has a species in the root's compartment, and whose substrate and product sets (H+ ignored) are not
those of a reaction already in that compartment. The model is only changed in memory.

Why two solves per candidate are enough. The candidate is added as a reversible reaction and its smallest and
largest flux are computed.
  - Both zero: the candidate must stay at zero, so every other flux range is what it was without it.
  - Otherwise: a root is a gap metabolite (all of its reactions are blocked). The candidate has a non-zero
    coefficient on the root, so the root's balance forces another of the root's reactions to carry flux, and
    that reaction was blocked before.
"""
import collections
import hashlib
import itertools
import json
import math
import os
import sys

import pandas as pd

H_PLUS = "MNXM1"
SOLVER_TOLERANCE = {"glpk": 1e-7, "gurobi": 1e-9}   # cobra's default zero cutoff; the census Gurobi setting
SPOT_CHECK_FAILS = 20                               # failing candidates re-counted as an implementation check

# Plan-time counts on 7ccbb568 for the 364 census roots (TASK.md). The candidate list must reproduce them.
EXPECTED_PAIRS, EXPECTED_PLACEMENTS = 391, 400
EXPECTED_FUNNEL = {
    "in no open reaction": 28,
    "no MetaNetX id": 61,
    "no MetaNetX enzyme reaction uses the compound": 81,
    "no reaction marked balanced": 35,
    "every balanced reaction needs a compound the compartment lacks": 58,
    "only candidates already in the compartment": 20,
    "one candidate": 23,
    "two or more candidates": 58,
}


def check_input_hashes(who, artifacts=""):
    """Stop unless every input listed for this script in input_hashes.csv has the recorded SHA-256.
    Returns {key: path}."""
    here = os.path.dirname(os.path.abspath(__file__))
    bases = {"repo": os.path.abspath(os.path.join(here, "..", "..", "..", "..")),
             "research": os.environ.get("IYALI26_RESEARCH_ROOT", ""), "artifacts": artifacts}
    table = pd.read_csv(os.path.join(here, "..", "input_hashes.csv"))
    paths = {}
    for row in table[table.checked_by.isin([who, "both"])].itertuples():
        if not bases[row.base]:
            sys.exit(f"GATE FAILED: no {row.base} folder to locate {row.key}")
        path = os.path.join(bases[row.base], row.path)
        with open(path, "rb") as handle:
            digest = hashlib.sha256(handle.read()).hexdigest()
        if digest != row.sha256:
            sys.exit(f"GATE FAILED: {row.key} has SHA-256 {digest}, expected {row.sha256}")
        paths[row.key] = path
    print(f"gate passed: {len(paths)} input hashes equal input_hashes.csv")
    return paths


# ------------------------------------------------------------------------------------------ candidate list
def load_universe(path):
    """MetaNetX non-transport reactions: id -> (net stoichiometry by compound id, balance mark)."""
    universe = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            rid, equation, _ref, _cls, balanced, transport = (line.rstrip("\n").split("\t") + [""] * 6)[:6]
            if rid == "EMPTY" or transport == "T":
                continue
            left, right = equation.split(" = ") if " = " in equation else (equation.replace(" =", ""), "")
            stoich = collections.defaultdict(float)
            for side, sign in ((left, -1), (right, 1)):
                for term in side.split(" + "):
                    term = term.strip()
                    if term:
                        coefficient, compound = term.split(" ", 1)
                        stoich[compound.split("@")[0]] += sign * float(coefficient)
            universe[rid] = ({k: v for k, v in stoich.items() if v != 0}, balanced == "B")
    return universe


def load_pool(gapfill_dir):
    """MetaNetX reaction id -> gene ids recorded in the gene-based gap-fill tables."""
    pool = collections.defaultdict(set)
    for name in ("gap_fill_candidates.csv", "gap_fill_prioritized.csv"):
        table = pd.read_csv(os.path.join(gapfill_dir, name), dtype=str).fillna("")
        for row in table.itertuples():
            if row.mnxr_id:
                pool[row.mnxr_id].add(row.gene_id)
    return pool


def signature(stoich):
    """Substrate set and product set with H+ ignored, without direction."""
    return frozenset((frozenset(k for k, v in stoich.items() if v < 0 and k != H_PLUS),
                      frozenset(k for k, v in stoich.items() if v > 0 and k != H_PLUS)))


def mnx_id(metabolite):
    value = metabolite.annotation.get("metanetx.chemical")
    if isinstance(value, (list, tuple)):
        if len(value) != 1:
            sys.exit(f"{metabolite.id} carries {len(value)} MetaNetX ids; one is expected")
        value = value[0]
    return value or None


def enumerate_candidates(model, roots, universe, pool):
    """Placements to test, pairs identical to an existing reaction, and one funnel row per root."""
    by_compound = collections.defaultdict(set)
    for rid, (stoich, _balanced) in universe.items():
        for compound in stoich:
            by_compound[compound].add(rid)
    ident = {met.id: mnx_id(met) for met in model.metabolites}
    species_of = collections.defaultdict(list)           # (compartment, MetaNetX id) -> species ids
    for met in model.metabolites:
        if ident[met.id]:
            species_of[(met.compartment, ident[met.id])].append(met.id)
    existing = collections.defaultdict(list)             # signature -> [(reaction id, compartment, closed)]
    for rxn in model.reactions:
        compartments = {met.compartment for met in rxn.metabolites}
        if len(compartments) != 1 or len(rxn.metabolites) < 2 or any(ident[met.id] is None for met in rxn.metabolites):
            continue
        stoich = collections.defaultdict(float)
        for met, coefficient in rxn.metabolites.items():
            stoich[ident[met.id]] += coefficient
        existing[signature({k: v for k, v in stoich.items() if v != 0})].append(
            (rxn.id, compartments.pop(), rxn.lower_bound == 0 and rxn.upper_bound == 0))

    placements, identical, funnel = [], [], []
    for root in roots.itertuples():
        compound, compartment = ident[root.id], root.compartment
        n_all = n_balanced = n_fit_balanced = 0
        kept = []
        if root.root_class.startswith("isolated"):
            status = "in no open reaction"
        elif compound is None:
            status = "no MetaNetX id"
        else:
            reactions = sorted(by_compound.get(compound, ()))
            n_all = len(reactions)
            n_balanced = sum(universe[rid][1] for rid in reactions)
            for rid in reactions:
                stoich, is_balanced = universe[rid]
                if any((compartment, c) not in species_of for c in stoich):
                    continue
                n_fit_balanced += is_balanced
                same = existing.get(signature(stoich), [])
                here = [s for s in same if s[1] == compartment]
                if here:
                    identical.append(dict(root=root.id, mnxr=rid, balanced=is_balanced,
                                          model_reactions=";".join(s[0] for s in here),
                                          any_open=any(not s[2] for s in here), pooled=root.pooled))
                    continue
                kept.append(rid)
                choices = [[root.id] if c == compound else species_of[(compartment, c)] for c in stoich]
                for variant, combo in enumerate(itertools.product(*choices)):
                    placements.append(dict(
                        root=root.id, mnxr=rid, variant=variant, balanced=is_balanced, pooled=root.pooled,
                        pool_genes=";".join(sorted(pool.get(rid, ()))),
                        copy_of=";".join(f"{s[0]}[{s[1]}]" for s in same),
                        root_coefficient=stoich[compound],
                        stoich=dict(zip(combo, (stoich[c] for c in stoich)))))
            kept_balanced = [rid for rid in kept if universe[rid][1]]
            if not reactions:
                status = "no MetaNetX enzyme reaction uses the compound"
            elif not n_balanced:
                status = "no reaction marked balanced"
            elif not n_fit_balanced:
                status = "every balanced reaction needs a compound the compartment lacks"
            elif not kept_balanced:
                status = "only candidates already in the compartment"
            elif len(kept_balanced) == 1:
                status = "one candidate"
            else:
                status = "two or more candidates"
        funnel.append(dict(root=root.id, name=root.name, compartment=compartment, root_class=root.root_class, um=root.um,
                           pooled=root.pooled, mnx_id=compound or "", status=status, n_reactions_with_compound=n_all,
                           n_marked_balanced=n_balanced, n_candidates=len(kept),
                           n_candidates_balanced=sum(universe[rid][1] for rid in kept)))
    return placements, identical, funnel


# ------------------------------------------------------------------------------------------------ LP tests
def load_model(path, solver):
    import cobra
    model = cobra.io.read_sbml_model(path)
    model.solver = solver
    if solver == "glpk":
        model.solver.configuration.presolve = True
    else:
        model.tolerance = SOLVER_TOLERANCE[solver]
        model.solver.configuration.verbosity = 0
        try:
            model.solver.problem.Params.Threads = 1
        except Exception:
            pass
    for rxn in model.exchanges:                          # the census "open" condition
        rxn.bounds = (-1000, 1000)
    return model


def add_candidate(model, stoich):
    import cobra
    rxn = cobra.Reaction("CANDIDATE_TEST", lower_bound=-1000, upper_bound=1000)
    model.add_reactions([rxn])
    rxn.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoich.items()})
    return rxn


def flux_range(model, reaction, counter):
    """Smallest and largest flux of one reaction, each with its raw solver status."""
    out = {}
    with model:
        model.objective = reaction
        for direction in ("min", "max"):
            model.objective_direction = direction
            value = model.slim_optimize(error_value=float("nan"))
            status = model.solver.status
            counter["lps"] += 1
            if status != "optimal" and model.solver.interface.__name__.endswith("glpk_interface"):
                model.solver.configuration.presolve = False       # the census scripts' single retry
                value = model.slim_optimize(error_value=float("nan"))
                status = model.solver.status + " (retried without presolve)"
                model.solver.configuration.presolve = True
                counter["lps"] += 1
            out[direction] = (value, status)
    return out["min"][0], out["max"][0], out["min"][1], out["max"][1]


def works(low, high, status_low, status_high, tolerance):
    """True or False; None when a solve did not end optimal (unresolved, never counted as a failure)."""
    if not (status_low.startswith("optimal") and status_high.startswith("optimal")) or math.isnan(low) or math.isnan(high):
        return None
    return max(abs(low), abs(high)) > tolerance


def newly_active(model, reaction_ids, counter):
    """Reactions of the list that are no longer blocked (FVA), or None when the solver fails twice."""
    from cobra.flux_analysis import find_blocked_reactions
    reactions = [model.reactions.get_by_id(r) for r in reaction_ids]
    kwargs = dict(reaction_list=reactions, open_exchanges=False, processes=1)
    counter["lps"] += 2 * len(reactions)
    try:
        still = set(find_blocked_reactions(model, **kwargs))
    except Exception:
        model.solver.configuration.presolve = False
        try:
            counter["lps"] += 2 * len(reactions)
            still = set(find_blocked_reactions(model, **kwargs))
        except Exception:
            return None
        finally:
            model.solver.configuration.presolve = True
    return [r for r in reaction_ids if r not in still]


def run_solver(model_path, solver, placements, census_roots, module_reactions, root_module, root_reactions):
    """Both passes for one solver. Returns one result dict per placement and the LP count."""
    model = load_model(model_path, solver)
    tolerance = SOLVER_TOLERANCE[solver]
    counter = dict(lps=0)
    results = [dict() for _ in placements]

    spot_left = SPOT_CHECK_FAILS if solver == "glpk" else 0
    for p, res in zip(placements, results):                       # pass 1: the candidate alone
        with model:
            rxn = add_candidate(model, p["stoich"])
            low, high, s_low, s_high = flux_range(model, rxn, counter)
            ok = works(low, high, s_low, s_high, tolerance)
            res.update({f"{solver}_alone_min": low, f"{solver}_alone_max": high,
                        f"{solver}_alone_status": f"{s_low} / {s_high}", f"works_alone_{solver}": ok})
            if solver != "glpk":
                continue
            module = module_reactions.get(root_module[p["root"]], [])
            res.update(n_module_blocked=len(module), n_unblocked_alone=None, unblocked_alone="",
                       unblocked_touches_root=None, spot_check_unblocked=None)
            if ok:
                unblocked = newly_active(model, module, counter)
                res["n_unblocked_alone"] = None if unblocked is None else len(unblocked)
                res["unblocked_alone"] = "" if unblocked is None else ";".join(unblocked)
                res["unblocked_touches_root"] = None if unblocked is None else bool(set(unblocked) & root_reactions[p["root"]])
            elif ok is False and spot_left and module and len(module) <= 50:
                unblocked = newly_active(model, module, counter)
                res["spot_check_unblocked"] = None if unblocked is None else len(unblocked)
                spot_left -= 1

    with model:                                                   # pass 2: every other census root opened
        import cobra
        sinks = {}
        for i, root_id in enumerate(census_roots):
            sink = cobra.Reaction(f"OPEN_ROOT_{i}", lower_bound=-1000, upper_bound=1000)
            model.add_reactions([sink])
            sink.add_metabolites({model.metabolites.get_by_id(root_id): -1})
            sinks[root_id] = sink
        for p, res in zip(placements, results):
            with model:
                if p["root"] in sinks:
                    sinks[p["root"]].bounds = (0, 0)
                rxn = add_candidate(model, p["stoich"])
                low, high, s_low, s_high = flux_range(model, rxn, counter)
                res.update({f"{solver}_open_min": low, f"{solver}_open_max": high,
                            f"{solver}_open_status": f"{s_low} / {s_high}",
                            f"works_open_{solver}": works(low, high, s_low, s_high, tolerance)})
    return results, counter["lps"]


def direction(low, high, tolerance):
    forward, backward = high > tolerance, low < -tolerance
    return "both" if forward and backward else "forward" if forward else "backward" if backward else "none"


def stage_test(model_path, art, reac_prop, gapfill_dir):
    recorded = check_input_hashes("candidate_tests")
    given = {"model": model_path, "reac_prop": reac_prop,
             "gap_fill_candidates": os.path.join(gapfill_dir, "gap_fill_candidates.csv"),
             "gap_fill_prioritized": os.path.join(gapfill_dir, "gap_fill_prioritized.csv")}
    for key, path in given.items():
        if not os.path.samefile(path, recorded[key]):
            sys.exit(f"GATE FAILED: the {key} argument {path} is not the file listed in input_hashes.csv")
    out, tag = os.path.join(art, "out"), "7ccbb568"
    reactions = pd.read_csv(f"{out}/{tag}_open_reactions.csv")
    reactions["um"] = reactions["um"].fillna("")
    mets = pd.read_csv(f"{out}/{tag}_open_metabolites.csv").fillna("")
    met_info = mets.set_index("id")
    roots = pd.read_csv(f"{out}/{tag}_open_roots.csv").fillna("")[["id", "name", "compartment", "root_class", "um"]]
    roots["pooled"] = False
    pooled = pd.read_csv(f"{out}/{tag}_open_pooled_roots.csv")
    extra = []
    for row in pooled[~pooled.already_root].itertuples():
        for copy in str(row.copies).split(";"):
            extra.append(dict(id=copy, name=met_info.at[copy, "name"], compartment=met_info.at[copy, "compartment"],
                              root_class=row.pooled_class, um=met_info.at[copy, "um"], pooled=True))
    all_roots = pd.concat([roots, pd.DataFrame(extra)], ignore_index=True).sort_values("id", kind="stable")
    not_gap = [r for r in all_roots.id if not bool(met_info.at[r, "gap"])]
    if not_gap:
        sys.exit(f"GATE FAILED: roots that are not gap metabolites: {not_gap}")

    universe, pool = load_universe(reac_prop), load_pool(gapfill_dir)
    placements, identical, funnel = enumerate_candidates(load_model(model_path, "glpk"), all_roots, universe, pool)
    funnel = pd.DataFrame(funnel)
    census = [p for p in placements if not p["pooled"]]
    got = dict(pairs=len({(p["root"], p["mnxr"]) for p in census}), placements=len(census),
               funnel=funnel[~funnel.pooled].status.value_counts().to_dict())
    if got != dict(pairs=EXPECTED_PAIRS, placements=EXPECTED_PLACEMENTS, funnel=EXPECTED_FUNNEL):
        sys.exit(f"GATE FAILED: candidate list {got} differs from the plan-time count")
    print("gate passed: candidate list equals the plan-time count;", len(placements) - len(census), "placements at pooled roots")

    module_reactions = reactions[reactions.blocked & (reactions.um != "")].groupby("um").id.apply(list).to_dict()
    root_module = dict(zip(all_roots.id, all_roots.um))
    root_reactions = {r: set(str(met_info.at[r, "reactions"]).split(";")) for r in all_roots.id}
    census_roots = list(roots.id)
    merged = [dict(p) for p in placements]
    lps = {}
    for solver in ("glpk", "gurobi"):
        results, lps[solver] = run_solver(model_path, solver, placements, census_roots, module_reactions,
                                          root_module, root_reactions)
        for row, res in zip(merged, results):
            row.update(res)
        print(solver, "done:", lps[solver], "LPs", flush=True)

    names = met_info.name.to_dict()
    for row in merged:
        stoich = row.pop("stoich")
        left = " + ".join(f"{-c:g} {names[m]}" for m, c in stoich.items() if c < 0)
        right = " + ".join(f"{c:g} {names[m]}" for m, c in stoich.items() if c > 0)
        row["equation"] = f"{left} <=> {right}"
        row["species"] = json.dumps(stoich, sort_keys=True)
        for which in ("alone", "open"):
            row[f"direction_{which}"] = direction(row[f"glpk_{which}_min"], row[f"glpk_{which}_max"], SOLVER_TOLERANCE["glpk"])
    tests = pd.DataFrame(merged)

    spot = tests["spot_check_unblocked"].dropna().astype(float)
    passed = tests[tests.works_alone_glpk == True]  # noqa: E712
    problems = []
    if (spot > 0).any():
        problems.append(f"{int((spot > 0).sum())} failing candidates unblocked a reaction in the spot check")
    if (passed.n_unblocked_alone == 0).any() or (passed.unblocked_touches_root == False).any():  # noqa: E712
        problems.append("a working candidate left every reaction of the root blocked")
    if passed.n_unblocked_alone.isna().any():
        problems.append(f"{int(passed.n_unblocked_alone.isna().sum())} re-counts for working candidates did not solve")
    tests.to_csv(f"{out}/candidate_tests.csv", index=False)
    funnel.to_csv(f"{out}/candidate_funnel.csv", index=False)
    pd.DataFrame(identical, columns=["root", "mnxr", "balanced", "model_reactions", "any_open", "pooled"]).to_csv(
        f"{out}/candidate_identical.csv", index=False)
    summary = dict(placements=len(tests), census_placements=len(census), pairs=got["pairs"], lps=lps,
                   spot_checked_fails=int(len(spot)), spot_check_unblocked_max=float(spot.max()) if len(spot) else 0.0,
                   consistency_problems=problems)
    json.dump(summary, open(f"{out}/candidate_summary.json", "w"), indent=1)
    print(json.dumps(summary))
    if problems:
        sys.exit("CONSISTENCY CHECK FAILED: " + "; ".join(problems))


# ---------------------------------------------------------------------------------------------- tables
def stage_tables(art, report, decisions_path):
    out, tag = os.path.join(art, "out"), "7ccbb568"
    tests = pd.read_csv(f"{out}/candidate_tests.csv").fillna({"pool_genes": "", "copy_of": "", "unblocked_alone": ""})
    funnel = pd.read_csv(f"{out}/candidate_funnel.csv").fillna("")
    identical = pd.read_csv(f"{out}/candidate_identical.csv")
    mets = pd.read_csv(f"{out}/{tag}_open_metabolites.csv").fillna("").set_index("id")
    fix = json.load(open(f"{out}/{tag}_open_fixtests_out.json"))
    decisions = pd.read_csv(decisions_path)
    problems = json.load(open(f"{out}/candidate_summary.json"))["consistency_problems"]
    if problems:
        sys.exit("GATE FAILED: the test stage recorded: " + "; ".join(problems))
    failed = [t for t in fix if "solver failed" in t.get("note", "")]
    if failed:
        sys.exit(f"GATE FAILED: {len(failed)} census fix tests ended with a solver failure, first: "
                 f"{failed[0]['kind']} {failed[0]['root']}")

    def fix_result(kind):
        return {t["root"]: f"{t['n_unblocked']}/{t['n_um_rxns']}" for t in fix if t["kind"] == kind}
    sink, transport, exchange = (fix_result(k) for k in ("sink (source/sink) on root", "add transport", "add exchange"))

    def label(root):
        mine = set(str(mets.at[root, "reactions"]).split(";"))
        return "; ".join(d.decision for d in decisions.itertuples() if mine & set(d.reaction_ids.split(";")))

    def count(frame, column):
        """Distinct candidate reactions of one root for which any placement is True in the column."""
        return frame[frame[column] == True].mnxr.nunique()  # noqa: E712

    rows = []
    for f in funnel.itertuples():
        mine = tests[tests.root == f.root]
        bal = mine[mine.balanced]
        resolved = ~mine[["works_alone_glpk", "works_alone_gurobi", "works_open_glpk", "works_open_gurobi"]].isna().any(axis=1)
        unresolved = int((~resolved).sum())
        disagree = int((resolved & ((mine.works_alone_glpk != mine.works_alone_gurobi)
                                    | (mine.works_open_glpk != mine.works_open_gurobi))).sum())
        rows.append(dict(
            root=f.root, name=f.name, compartment=f.compartment, root_class=f.root_class, um=f.um, pooled=f.pooled,
            mnx_id=f.mnx_id, status=f.status, n_candidates=f.n_candidates, n_candidates_balanced=f.n_candidates_balanced,
            works_alone=count(mine, "works_alone_glpk"), works_alone_balanced=count(bal, "works_alone_glpk"),
            works_open=count(mine, "works_open_glpk"), works_open_balanced=count(bal, "works_open_glpk"),
            works_open_balanced_gurobi=count(bal, "works_open_gurobi"),
            works_open_in_gene_pool=count(mine[mine.pool_genes != ""], "works_open_glpk"),
            most_unblocked_alone=int(mine.n_unblocked_alone.max()) if mine.n_unblocked_alone.notna().any() else 0,
            placements_solvers_disagree=disagree, placements_unresolved=unresolved,
            sink_test=sink.get(f.root, ""), transport_test=transport.get(f.root, ""), exchange_test=exchange.get(f.root, ""),
            recorded_decision=label(f.root)))
    per_root = pd.DataFrame(rows)

    def bins(series):
        return {"none": int((series == 0).sum()), "one": int((series == 1).sum()), "two or more": int((series >= 2).sum())}
    census = per_root[~per_root.pooled]
    with_candidates = census[census.n_candidates > 0]
    # A root with an unresolved solve is not binned as "none works": it gets its own column.
    resolved = census[census.placements_unresolved == 0]
    summary = []
    for label_, column in (("candidates listed (balance-marked)", "n_candidates_balanced"), ("candidates listed (all)", "n_candidates"),
                           ("work alone (balance-marked)", "works_alone_balanced"), ("work alone (all)", "works_alone"),
                           ("work with other roots opened (balance-marked)", "works_open_balanced"),
                           ("work with other roots opened (all)", "works_open"),
                           ("work with other roots opened (balance-marked, Gurobi)", "works_open_balanced_gurobi")):
        frame = census if column.startswith("n_candidates") else resolved
        summary.append(dict(measure=label_, roots_total=len(census), **bins(frame[column]),
                            unresolved=len(census) - len(frame)))
    summary = pd.DataFrame(summary)

    test_columns = ["root", "mnxr", "variant", "balanced", "pooled", "pool_genes", "copy_of", "root_coefficient", "equation",
                    "works_alone_glpk", "works_alone_gurobi", "works_open_glpk", "works_open_gurobi", "direction_alone",
                    "direction_open", "n_module_blocked", "n_unblocked_alone", "glpk_alone_min", "glpk_alone_max",
                    "gurobi_alone_min", "gurobi_alone_max", "glpk_open_min", "glpk_open_max", "gurobi_open_min",
                    "gurobi_open_max", "glpk_alone_status", "gurobi_alone_status", "glpk_open_status", "gurobi_open_status"]
    tests[test_columns].to_csv(f"{report}/candidate_tests.csv", index=False)
    per_root.to_csv(f"{report}/root_candidate_counts.csv", index=False)
    identical.to_csv(f"{report}/candidates_identical_to_model.csv", index=False)
    summary.to_csv(f"{report}/candidate_summary.csv", index=False)
    os.makedirs(f"{art}/tables", exist_ok=True)
    tests.to_csv(f"{art}/tables/candidate_tests_full.csv", index=False)

    decided = census[census.recorded_decision != ""]
    print(json.dumps(dict(
        roots=len(census), pooled_copies=int(per_root.pooled.sum()),
        status=census.status.value_counts().to_dict(),
        placements=len(tests), placements_census=int((~tests.pooled).sum()),
        works_alone_placements=int((tests.works_alone_glpk == True).sum()),  # noqa: E712
        works_open_placements=int((tests.works_open_glpk == True).sum()),  # noqa: E712
        solver_disagreements=int(per_root.placements_solvers_disagree.sum()), unresolved=int(per_root.placements_unresolved.sum()),
        summary=summary.set_index("measure")[["none", "one", "two or more", "unresolved"]].to_dict("index"),
        reading_count_roots_two_or_more_open_balanced=int((resolved.works_open_balanced >= 2).sum()),
        roots_with_an_unresolved_solve=len(census) - len(resolved),
        roots_with_candidates=len(with_candidates),
        roots_under_recorded_decision=len(decided),
        roots_under_recorded_decision_with_working_candidate=int((decided.works_open > 0).sum())), indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        stage_test(*sys.argv[2:6])
    elif len(sys.argv) > 1 and sys.argv[1] == "tables":
        stage_tables(*sys.argv[2:5])
    else:
        sys.exit(__doc__)
