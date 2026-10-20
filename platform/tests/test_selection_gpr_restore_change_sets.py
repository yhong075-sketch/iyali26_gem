"""In-memory checks for the R671 and CPA/URA2 selection restores and the corrected
mis-annotation removal (2026-10-19); no optimization."""

import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.model_layout import MODEL, PLATFORM_ROOT
from scripts.gem_annotate.patches import remove_misannotated_gprs
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.selection_gpr_restore import (
    CHANGE_SETS, CLEARED_GPR_RESTORE, CPA_URA2_RESTORE, R671_GPR_RESTORE, TARGETS,
    apply_selection_gpr_restore, restore_after_selection)

BASELINE = "E5_coq9_resp_pkg_r132_20261012.xml"
NEW_CHANGE_SETS = (R671_GPR_RESTORE, CPA_URA2_RESTORE)
SPECS = {change_id: json.loads(CHANGE_SETS[change_id][0].read_text()) for change_id in NEW_CHANGE_SETS}
# Applies every change set to its cleared fixture and prints each target's species order.
ORDER_PROBE = """
import json, sys
sys.path.insert(0, "tests")
import test_selection_gpr_restore as cleared, test_selection_gpr_restore_change_sets as sets
orders = {}
for change_id, model in [(sets.CLEARED_GPR_RESTORE, cleared.fixture())] + [(c, sets.fixture(c)) for c in sets.NEW_CHANGE_SETS]:
    sets.apply_selection_gpr_restore(model, change_id=change_id)
    orders.update({rid: [m.id for m in model.reactions.get_by_id(rid).metabolites] for rid in sets.CHANGE_SETS[change_id][1]})
print(json.dumps(orders))
"""


def _add(model, rid, state):
    reaction = Reaction(rid, lower_bound=state["bounds"][0], upper_bound=state["bounds"][1])
    model.add_reactions([reaction])
    for mid in state["stoichiometry"]:
        if mid not in model.metabolites:
            model.add_metabolites([Metabolite(mid, compartment="C_cy")])
    reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in state["stoichiometry"].items()})
    reaction.gene_reaction_rule = state["gpr"]
    return reaction


def fixture(change_id, state="cleared"):
    spec = SPECS[change_id]
    model = Model("selection_restore_change_set_fixture")
    model.add_metabolites([Metabolite(mid, formula=row["formula"], charge=row["charge"],
                                      compartment=row["compartment"])
                           for mid, row in spec["species"].items()])
    context = spec.get("context", {})
    for rid, expected in context.items():
        if rid not in spec["reactions"]:
            _add(model, rid, expected)
    for rid, row in spec["reactions"].items():
        reaction = _add(model, rid, row[state])
        reaction.name = context.get(rid, {}).get("name", "")
        if "ec-code" in context.get(rid, {}):
            reaction.annotation["ec-code"] = list(context[rid]["ec-code"])
    for row in spec["reactions"].values():
        for gene in row["genes"]:
            model.genes.get_by_id(gene["id"]).annotation.update(
                {k: v for k, v in gene.items() if k not in ("id", "name")})
    return model


def snapshot(model):
    return model_definition(model), {r.id: dict(r.notes) for r in model.reactions}


class ChangeSetRestoreTests(unittest.TestCase):
    def test_restores_targets_and_is_idempotent(self):
        for change_id, spec in SPECS.items():
            with self.subTest(change_id):
                model = fixture(change_id)
                genes, metabolites = {g.id for g in model.genes}, {m.id for m in model.metabolites}
                result = apply_selection_gpr_restore(model, change_id=change_id)
                self.assertEqual(result["status"], "applied")
                for rid, row in spec["reactions"].items():
                    reaction = model.reactions.get_by_id(rid)
                    self.assertEqual({m.id: c for m, c in reaction.metabolites.items()},
                                     row["restored"]["stoichiometry"])
                    self.assertEqual(reaction.gene_reaction_rule, row["restored"]["gpr"])
                    self.assertEqual(reaction.check_mass_balance(), row.get("unchanged_equation_residual", {}))
                    self.assertEqual(reaction.notes["selection_gpr_restore"], change_id)
                self.assertEqual({g.id for g in model.genes}, genes)  # No gene added.
                self.assertEqual({m.id for m in model.metabolites}, metabolites)  # No metabolite added.
                state = snapshot(model)
                self.assertEqual(apply_selection_gpr_restore(model, change_id=change_id)["status"], "already_correct")
                self.assertEqual(snapshot(model), state)

    def test_r671_keeps_only_its_recorded_residual(self):
        result = apply_selection_gpr_restore(fixture(R671_GPR_RESTORE), change_id=R671_GPR_RESTORE)
        self.assertEqual(result["reactions"][0]["unchanged_equation_residual"], {"charge": 1.0})

    def test_gpr_propagation(self):
        r671 = fixture(R671_GPR_RESTORE)
        apply_selection_gpr_restore(r671, change_id=R671_GPR_RESTORE)
        cpa = fixture(CPA_URA2_RESTORE)
        apply_selection_gpr_restore(cpa, change_id=CPA_URA2_RESTORE)
        cases = ((r671, "YALI1E11370g", {"R671": False}), (r671, "YALI1F23441g", {"R671": True}),
                 (cpa, "YALI1D09420g", {"R190": True, "R159": False}),
                 (cpa, "YALI1C33005g", {"R190": True, "R159": False}),
                 (cpa, "YALI1E11768g", {"R190": False, "R159": True}))
        for model, gene, closed in cases:
            with self.subTest(gene), model:
                model.genes.get_by_id(gene).knock_out()
                for rid, is_closed in closed.items():
                    self.assertEqual(model.reactions.get_by_id(rid).bounds == (0, 0), is_closed, rid)

    def test_r159_leaves_the_shared_carbamoyl_phosphate_pool(self):
        model = fixture(CPA_URA2_RESTORE)
        apply_selection_gpr_restore(model, change_id=CPA_URA2_RESTORE)
        pool = model.metabolites.get_by_id("m325[C_cy]")
        self.assertEqual(sorted(r.id for r in pool.reactions), ["R190", "R607"])

    def test_restored_species_order_does_not_depend_on_the_hash_seed(self):
        orders = []
        for seed in ("0", "1", "2"):
            completed = subprocess.run([sys.executable, "-B", "-c", ORDER_PROBE], cwd=PLATFORM_ROOT, text=True,
                                       capture_output=True, check=False, timeout=120,
                                       env={**os.environ, "PYTHONHASHSEED": seed})
            self.assertEqual(completed.returncode, 0, completed.stderr)
            orders.append(json.loads(completed.stdout.splitlines()[-1]))
        self.assertEqual(orders[1:], orders[:1] * 2)
        for change_id, (path, targets) in CHANGE_SETS.items():
            spec = json.loads(path.read_text())
            for rid in targets:
                cleared, restored = (spec["reactions"][rid][s]["stoichiometry"] for s in ("cleared", "restored"))
                # Kept species stay in place; added species follow in sorted ID order.
                expected = [m for m in cleared if m in restored] + sorted(set(restored) - set(cleared))
                self.assertEqual(orders[0][rid], expected, rid)

    def test_rerun_after_sbml_round_trip_is_a_no_op(self):
        for change_id in NEW_CHANGE_SETS:
            with self.subTest(change_id):
                model = fixture(change_id)
                apply_selection_gpr_restore(model, change_id=change_id)
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "restored.xml"
                    write_deterministic_sbml_model(model, path)
                    loaded = read_sbml_model(str(path))
                state = snapshot(loaded)
                self.assertEqual(apply_selection_gpr_restore(loaded, change_id=change_id)["status"], "already_correct")
                self.assertEqual(snapshot(loaded), state)

    def test_stale_inputs_fail_closed_without_edits(self):
        def r671_gpr(m): m.reactions.R671.gene_reaction_rule = "YALI1E11370g"
        def r671_equation(m): m.reactions.R671.add_metabolites({m.metabolites.get_by_id("m10[C_cy]"): -1.0})
        def nadp_charge(m): m.metabolites.get_by_id("m41[C_cy]").charge = 1
        def gene_identity(m): m.genes.YALI1F23441g.annotation["ncbigene"] = "0"
        def note(m): m.reactions.R671.notes["selection_gpr_restore"] = "other change"
        def r607_pool(m): m.reactions.R607.add_metabolites({m.metabolites.get_by_id("m325[C_cy]"): 1.0})
        def r607_gpr(m): m.reactions.R607.gene_reaction_rule = "YALI1E11768g"
        def r159_name(m): m.reactions.R159.name = "aspartate carbamoyltransferase"
        def r190_ec(m): m.reactions.R190.annotation["ec-code"] = ["6.3.5.5"]
        def mixed_state(m): m.reactions.R159.gene_reaction_rule = "YALI1E11768g"
        def cpa_identity(m): m.genes.YALI1D09420g.annotation["uniprot"] = "conflicting_accession"
        def missing(m): m.remove_reactions([m.reactions.R190])
        def extra_pool_user(m):
            m.add_reactions([Reaction("R_other")])
            m.reactions.R_other.add_metabolites({m.metabolites.get_by_id("m325[C_cy]"): -1.0})
        neither = "neither the cleared nor the restored state"
        cases = [(R671_GPR_RESTORE, r671_gpr, neither), (R671_GPR_RESTORE, r671_equation, neither),
                 (R671_GPR_RESTORE, nadp_charge, "species differs"),
                 (R671_GPR_RESTORE, gene_identity, "identity differs"),
                 (R671_GPR_RESTORE, note, "conflicting restore notes"),
                 (CPA_URA2_RESTORE, r607_pool, "R607: context stoichiometry"),
                 (CPA_URA2_RESTORE, r607_gpr, "R607: context gpr"),
                 (CPA_URA2_RESTORE, r159_name, "R159: context name"),
                 (CPA_URA2_RESTORE, r190_ec, "R190: context ec-code"),
                 (CPA_URA2_RESTORE, mixed_state, neither),
                 (CPA_URA2_RESTORE, cpa_identity, "identity differs"),
                 (CPA_URA2_RESTORE, missing, "R190: context reaction missing"),
                 (CPA_URA2_RESTORE, extra_pool_user, "reactions using the pool differ")]
        for change_id, change, message in cases:
            with self.subTest(change.__name__):
                model = fixture(change_id)
                change(model)
                before = snapshot(model)
                with self.assertRaisesRegex(ValueError, message):
                    apply_selection_gpr_restore(model, change_id=change_id)
                self.assertEqual(snapshot(model), before)

    def test_unexpected_records_are_rejected(self):
        no_residual = copy.deepcopy(SPECS[R671_GPR_RESTORE])
        del no_residual["reactions"]["R671"]["unchanged_equation_residual"]
        wrong_residual = copy.deepcopy(SPECS[R671_GPR_RESTORE])
        wrong_residual["reactions"]["R671"]["unchanged_equation_residual"] = {"charge": 2.0}
        edited_with_residual = copy.deepcopy(SPECS[CPA_URA2_RESTORE])
        edited_with_residual["reactions"]["R159"]["unchanged_equation_residual"] = {"H": 1.0}
        genes_mismatch = copy.deepcopy(SPECS[CPA_URA2_RESTORE])
        genes_mismatch["reactions"]["R190"]["genes"] = genes_mismatch["reactions"]["R190"]["genes"][:1]
        bad_id = copy.deepcopy(SPECS[CPA_URA2_RESTORE]); bad_id["change_id"] = CLEARED_GPR_RESTORE
        missing = copy.deepcopy(SPECS[CPA_URA2_RESTORE]); del missing["reactions"]["R190"]
        for name, change_id, spec, message in (
                ("no_residual", R671_GPR_RESTORE, no_residual, "not balanced"),
                ("wrong_residual", R671_GPR_RESTORE, wrong_residual, "not balanced"),
                ("edited_with_residual", CPA_URA2_RESTORE, edited_with_residual, "only an unchanged equation"),
                ("genes_mismatch", CPA_URA2_RESTORE, genes_mismatch, "GPR and genes disagree"),
                ("change_id", CPA_URA2_RESTORE, bad_id, "Unexpected selection GPR restore scope"),
                ("scope", CPA_URA2_RESTORE, missing, "Unexpected selection GPR restore scope")):
            with self.subTest(name):
                model = fixture(change_id)
                before = snapshot(model)
                with self.assertRaisesRegex(ValueError, message):
                    apply_selection_gpr_restore(model, spec, change_id=change_id)
                self.assertEqual(snapshot(model), before)

    def test_flags_are_opt_in_and_need_a_complete_selection(self):
        args = parse_args([])
        self.assertFalse(args.restore_selection_overwritten_gprs)
        self.assertFalse(args.restore_selection_cpa_ura2)
        args = parse_args(["--restore-selection-overwritten-gprs", "--restore-selection-cpa-ura2"])
        self.assertTrue(args.restore_selection_overwritten_gprs and args.restore_selection_cpa_ura2)
        for change_id in NEW_CHANGE_SETS:
            with self.subTest(change_id):
                model = fixture(change_id)
                before = snapshot(model)
                with self.assertRaises(ValueError):
                    restore_after_selection(model, {"complete": False}, change_id)
                self.assertEqual(snapshot(model), before)
                self.assertEqual(restore_after_selection(model, {"complete": True}, change_id)["status"], "applied")

    def test_baseline_candidate_changes_only_the_targets(self):
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        before = model_definition(model)
        for change_id in (CLEARED_GPR_RESTORE,) + NEW_CHANGE_SETS:
            self.assertEqual(apply_selection_gpr_restore(model, change_id=change_id)["status"], "applied")
        after = model_definition(model)
        self.assertEqual({k: v for k, v in after.items() if k != "reactions"},
                         {k: v for k, v in before.items() if k != "reactions"})
        changed = sorted(r for r in after["reactions"] if after["reactions"][r] != before["reactions"][r])
        self.assertEqual(changed, sorted(TARGETS + ("R671", "R159", "R190")))
        for rid in ("R159", "R190"):
            self.assertEqual(model.reactions.get_by_id(rid).check_mass_balance(), {})


class MisannotationRemovalTests(unittest.TestCase):
    def test_keeps_transketolase_and_removes_gatb(self):
        model = Model("misannotation_fixture")
        for rid, rule in (("R765", "YALI1E07744g or YALI1D02625g"), ("R766", "YALI1E07744g or YALI1D02625g"),
                          ("R671", "YALI1F23441g or YALI1E11370g")):
            reaction = Reaction(rid)
            model.add_reactions([reaction])
            reaction.gene_reaction_rule = rule
        self.assertEqual(remove_misannotated_gprs(model), 1)
        self.assertEqual(model.reactions.R765.gene_reaction_rule, "YALI1E07744g or YALI1D02625g")
        self.assertEqual(model.reactions.R766.gene_reaction_rule, "YALI1E07744g or YALI1D02625g")
        self.assertEqual(model.reactions.R671.gene_reaction_rule, "YALI1F23441g")
        self.assertEqual(remove_misannotated_gprs(model), 0)


if __name__ == "__main__":
    unittest.main()
