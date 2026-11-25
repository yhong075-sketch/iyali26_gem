"""Focused curation checks against a source-derived, eight-reaction SBML excerpt."""

import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra import Metabolite
from cobra.core.gene import GPR
from cobra.io import read_sbml_model, write_sbml_model
from optlang import symbolics
from swiglpk import (doubleArray, glp_get_col_name, glp_get_mat_row, glp_get_num_cols, glp_get_num_rows,
                     glp_get_row_name, intArray)

from scripts.gem_annotate.coq9 import (
    apply_coq9_curation, boolean_key, deduplicate_and, exact_residual, gene_evidence,
)
from scripts.gem_annotate.metabolites import normalize_all_annotations
from scripts.gem_annotate.patches import clean_ec_overload

FIXTURE = Path(__file__).parent / "fixtures" / "coq9_reference_extract.xml"


def expression_terms(expression):
    # One (variable name, coefficient, lower bound, upper bound, type) per term: what the "variables" entry of
    # semantics() holds for the variable, so an entry compared alone sees a change of its variables. The bounds are
    # the variable's numbers, so a bound written 0 in one model and 0.0 in another compares equal. A constant term
    # has none of the three.
    return sorted((term.name, float(coefficient), term.lb, term.ub, term.type) if hasattr(term, "name")
                  else (str(term), float(coefficient), None, None, None)
                  for term, coefficient in expression.as_coefficients_dict().items())


def constraint_terms(model):
    """Each constraint's terms, read from the solver's own matrix.

    Same content as expression_terms(constraint.expression) for every constraint, without building one symbolic
    expression per constraint (that took most of a snapshot's time).
    """
    model.solver.update()  # apply optlang's queued additions and removals, as reading model.constraints does
    problem, interface = model.solver.problem, model.problem.__name__
    rows = {}
    if interface == "optlang.gurobi_interface":
        for constraint in problem.getConstrs():
            row, name = problem.getRow(constraint), constraint.ConstrName
            terms = ((row.getVar(i).VarName, row.getCoeff(i)) for i in range(row.size()))
            # optlang gives a two-sided constraint a helper column "<name>_aux" and leaves it out of the expression.
            rows[name] = [term for term in terms if term[0] != name + "_aux"]
    elif interface == "optlang.glpk_interface":
        count = glp_get_num_cols(problem)
        columns = [None] + [glp_get_col_name(problem, j) for j in range(1, count + 1)]
        indices, values = intArray(count + 1), doubleArray(count + 1)
        for i in range(1, glp_get_num_rows(problem) + 1):
            length = glp_get_mat_row(problem, i, indices, values)
            rows[glp_get_row_name(problem, i)] = [(columns[indices[k]], values[k]) for k in range(1, length + 1)]
    else:
        return {c.name: expression_terms(c.expression) for c in model.constraints}
    variables = {v.name: (v.lb, v.ub, v.type) for v in model.variables}
    # A row without terms reads as the symbolic path writes it (the constant 0), so every path gives the same rows.
    empty = expression_terms(symbolics.add([]))
    return {name: sorted((variable, coefficient, *variables[variable]) for variable, coefficient in terms)
            or list(empty) for name, terms in rows.items()}


def semantics(model):
    """S, species, bounds/directions, objective, variables and ALL solver constraints."""
    rows = constraint_terms(model)
    return {
        "species": {m.id: (m.name, m.formula, m.charge, m.compartment) for m in model.metabolites},
        "compartments": model.compartments,
        "genes": sorted(g.id for g in model.genes),
        "reactions": {r.id: ({m.id: v for m, v in r.metabolites.items()}, r.bounds,
                               r.reversibility, boolean_key(r.gpr.body)) for r in model.reactions},
        "objective": (model.objective.direction, expression_terms(model.objective.expression)),
        "variables": {v.name: (v.lb, v.ub, v.type) for v in model.variables},
        "constraints": {c.name: (c.lb, c.ub, rows[c.name]) for c in model.constraints},
    }


def annotations(model):
    return {r.id: (r.name, copy.deepcopy(r.annotation), copy.deepcopy(r.notes), r.gene_reaction_rule)
            for r in model.reactions}


class CoQ9CurationTests(unittest.TestCase):
    def setUp(self):
        self.model = read_sbml_model(str(FIXTURE))

    def test_real_cli_default_and_all_modes_forward_to_same_builder(self):
        from scripts.gem_annotate.main import main
        with patch("scripts.gem_annotate.main.build_model") as build:
            main([])
            self.assertEqual(build.call_args.args[0].coq9_curation, "metadata")
            for mode in ("off", "metadata", "qcycle"):
                main(["--coq9-curation", mode, "--offline", "--no-solve"])
                self.assertEqual(build.call_args.args[0].coq9_curation, mode)
                self.assertTrue(build.call_args.args[0].no_solve)

    def test_metadata_preserves_full_mathematics_and_extra_constraint(self):
        m = self.model
        m.add_cons_vars(m.problem.Constraint(m.reactions.R305.flux_expression, ub=9, name="test_extra"))
        before = semantics(m)
        report = apply_coq9_curation(m)
        self.assertEqual(before, semantics(m))
        self.assertTrue(report["requested_mode_complete"])
        rows = {r.get("item"): r for r in report["items"]}
        self.assertEqual([rows[r]["status"] for r in ("R305", "R385", "R18", "R695", "R570", "R19", "R2062")], ["applied"] * 7)
        self.assertEqual((rows["R2062"]["before_occurrences"], rows["R2062"]["after_occurrences"]), (47, 28))

    def test_targets_history_and_export_readback(self):
        m = self.model
        original = semantics(m)
        unrelated = copy.deepcopy(m.reactions.R305.annotation["metanetx.reaction"])
        apply_coq9_curation(m)
        normalize_all_annotations(m)
        clean_ec_overload(m)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.xml"
            write_sbml_model(m, str(path))
            loaded = read_sbml_model(str(path))
        self.assertEqual(original, semantics(loaded))
        self.assertEqual(loaded.reactions.R305.name, "ubiquinol-9:cytochrome-c reductase (complex III)")
        for rid, ec in [("R305", "7.1.1.8"), ("R18", "2.1.1.201"), ("R695", "1.14.13.253"), ("R570", "1.6.5.9")]:
            self.assertIn(loaded.reactions.get_by_id(rid).annotation["ec-code"], (ec, [ec]))
            self.assertEqual(loaded.reactions.get_by_id(rid).notes["PROTEIN_CLASS"], ec)
        self.assertEqual(loaded.reactions.R385.notes["PROTEIN_CLASS"], "2.1.1.64")
        self.assertNotIn("ec-code", loaded.reactions.R19.annotation)
        self.assertNotIn("PROTEIN_CLASS", loaded.reactions.R19.notes)
        self.assertEqual(len(loaded.reactions.R19.notes["coq9_previous_ec_code"].split("; ")), 8)
        self.assertEqual(loaded.reactions.R305.annotation["metanetx.reaction"], unrelated)

    def test_wrong_r305_name_is_not_used_to_create_new_reaction_identity(self):
        from scripts.gem_annotate.reactions import annotate_reactions
        from scripts.gem_annotate.annotate_reactions_extended import annotate_remaining_reactions
        m = self.model
        r = m.reactions.R305
        old_name = r.name
        r.annotation = {}
        xrefs = {"by_mnxr": {"MNXR151780": [("ec-code", "1.9.3.1")]},
                 "desc_index": {old_name.lower(): "MNXR151780"}, "bigg_to_mnxr": {}}
        annotate_reactions(m, xrefs, name_exclusions={"R305": old_name})
        xrefs["ec_to_mnxr"] = {}
        annotate_remaining_reactions(m, xrefs, name_exclusions={"R305": old_name})
        self.assertNotIn("metanetx.reaction", r.annotation)
        report = apply_coq9_curation(m)
        self.assertTrue(report["requested_mode_complete"])
        self.assertEqual(r.annotation["ec-code"], ["7.1.1.8"])
        # A conflicting identity already in a different input is not erased.
        r.annotation["ec-code"] = ["1.9.3.1"]
        report = apply_coq9_curation(m)
        self.assertFalse(report["requested_mode_complete"])
        self.assertEqual(r.annotation["ec-code"], ["1.9.3.1"])

    def test_qcycle_only_two_coefficients_and_idempotency(self):
        m = self.model
        apply_coq9_curation(m)
        reference = m.copy()
        before = semantics(m)
        report = apply_coq9_curation(m, "qcycle")
        self.assertTrue(report["requested_mode_complete"])
        deltas = [(r.id, met.id, value - before["reactions"][r.id][0].get(met.id, 0))
                  for r in m.reactions for met, value in r.metabolites.items()
                  if value != before["reactions"][r.id][0].get(met.id, 0)]
        self.assertEqual(sorted(deltas), [("R305", "m10[C_cy]", 2.5), ("R305", "m28[C_mi]", -0.5)])
        self.assertEqual(exact_residual(m.reactions.R305), {})
        # Revert exactly those two entries in a separate copy; ALL constraints match.
        candidate = m.copy()
        candidate.reactions.R305.add_metabolites({candidate.metabolites.get_by_id("m10[C_cy]"): 1.5,
                                                candidate.metabolites.get_by_id("m28[C_mi]"): -1.5}, combine=False)
        self.assertEqual(semantics(candidate), semantics(reference))
        state = (semantics(m), annotations(m), gene_evidence(m))
        again = apply_coq9_curation(m, "qcycle")
        self.assertEqual(state, (semantics(m), annotations(m), gene_evidence(m)))
        self.assertEqual(next(r for r in again["items"] if r.get("item") == "R305_qcycle")["status"], "already_correct")
        self.assertIn("Q-cycle candidate present", m.reactions.R305.notes["coq9_proton_state"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "qcycle.xml"
            write_sbml_model(m, str(path))
            loaded = read_sbml_model(str(path))
        self.assertEqual(semantics(m), semantics(loaded))
        self.assertEqual(exact_residual(loaded.reactions.R305), {})
        readback_state = annotations(loaded)
        apply_coq9_curation(loaded, "qcycle")
        self.assertEqual(readback_state, annotations(loaded))

    def test_off_and_metadata_idempotency(self):
        m = self.model
        original = (semantics(m), annotations(m))
        apply_coq9_curation(m, "off")
        self.assertEqual(original, (semantics(m), annotations(m)))
        apply_coq9_curation(m)
        state = (semantics(m), annotations(m), gene_evidence(m))
        apply_coq9_curation(m)
        self.assertEqual(state, (semantics(m), annotations(m), gene_evidence(m)))

    def test_nested_boolean_deduplication_keeps_or_and_members(self):
        import ast
        rule = GPR.from_string("(a and (a and b)) or (c and c) or (d and (e or f) and (e or f))")
        reduced = GPR(ast.Expression(body=deduplicate_and(rule.body)))
        self.assertEqual(rule.genes, reduced.genes)
        self.assertEqual(boolean_key(rule.body), boolean_key(reduced.body))
        # Small independent Boolean check, never a 2**28 enumeration.
        for knockout in [set(), {"a"}, {"b", "c"}, {"a", "c", "e", "f"}, set(rule.genes)]:
            self.assertEqual(rule.eval(knockout), reduced.eval(knockout))
        repeated_or = GPR.from_string("a or a")
        self.assertEqual(ast.dump(repeated_or.body), ast.dump(deduplicate_and(repeated_or.body)))

    def test_local_conflicts_preserve_changed_reactions_and_other_items_continue(self):
        for mutation in ("species", "equation", "gpr", "bounds", "annotation"):
            with self.subTest(mutation=mutation):
                m = self.model.copy()
                r = m.reactions.R305
                if mutation == "species":
                    m.metabolites.get_by_id("m28[C_mi]").compartment = "C_cy"
                elif mutation == "equation":
                    r.add_metabolites({m.metabolites.get_by_id("m28[C_mi]"): -3}, combine=False)
                elif mutation == "gpr":
                    r.gene_reaction_rule += " or novel_gene"
                elif mutation == "bounds":
                    r.upper_bound = 7
                else:
                    r.annotation["ec-code"] = ["1.1.1.1"]
                before = (semantics(m)["reactions"]["R305"], annotations(m)["R305"])
                report = apply_coq9_curation(m, "qcycle")
                self.assertFalse(report["requested_mode_complete"])
                self.assertEqual(before, (semantics(m)["reactions"]["R305"], annotations(m)["R305"]))
                self.assertEqual(m.reactions.R385.notes["PROTEIN_CLASS"], "2.1.1.64")

    def test_evidence_and_native_overlay_are_observed_not_forced(self):
        m = self.model
        evidence = {r["yali1_gene"]: r for r in gene_evidence(m)}
        self.assertEqual(evidence["YALI1A08781g"]["representation_status"], "not represented / not testable")
        m.reactions.R385.gene_reaction_rule += " or YALI1A08781g"
        m.reactions.R573.upper_bound = 1000
        report = apply_coq9_curation(m)
        evidence = {r["yali1_gene"]: r for r in gene_evidence(m)}
        self.assertIn("R385", evidence["YALI1A08781g"]["actual_reactions"])
        self.assertTrue(evidence["YALI1A08781g"]["association_changed_from_reference"])
        self.assertEqual(m.reactions.R573.bounds, (0, 1000))
        self.assertFalse(report["requested_mode_complete"])


class SnapshotTests(unittest.TestCase):
    """semantics() reads the constraint rows from the solver's matrix: one reader for Gurobi and one for GLPK."""

    def model(self, solver):
        m = read_sbml_model(str(FIXTURE))
        m.solver = solver
        flux = m.reactions.R305.flux_expression
        m.add_cons_vars([m.problem.Constraint(flux, ub=9, name="one_sided"),
                         m.problem.Constraint(flux, lb=-3, ub=9, name="two_sided")])
        m.add_metabolites([Metabolite("unused[C_cy]", compartment="C_cy")])  # its row has no terms
        return m

    def matrix_only(self, m):
        """Fail when a constraint's symbolic expression is built: that is the path for other solvers, not a reader."""
        return patch.object(m.problem.Constraint, "expression",
                            property(lambda constraint: self.fail("a constraint row was built symbolically")))

    @staticmethod
    def changed(m, edit):
        """The entries of the snapshot, and the constraint rows, that differ after the edit."""
        before = semantics(m)
        edit()
        after = semantics(m)
        return ({key for key in after if after[key] != before[key]},
                {name for name, row in after["constraints"].items() if row != before["constraints"][name]})

    def test_matrix_rows_equal_the_symbolic_rows(self):
        for solver in ("gurobi", "glpk"):
            with self.subTest(solver=solver):
                m = self.model(solver)
                direct = constraint_terms(m)  # read first: optlang has not yet applied the last added row
                with self.matrix_only(m):
                    snapshot = {name: row[2] for name, row in semantics(m)["constraints"].items()}
                symbolic = {c.name: expression_terms(c.expression) for c in m.constraints}
                self.assertEqual((direct, snapshot), (symbolic, symbolic))
                r305 = m.reactions.R305
                # Each term: variable, coefficient, lower bound, upper bound, type.
                flux = [(r305.id, 1.0, 0, 1000, "continuous"), (r305.reverse_id, -1.0, 0, 0, "continuous")]
                self.assertEqual((symbolic["one_sided"], symbolic["two_sided"]), (flux, flux))
                self.assertEqual(symbolic["unused[C_cy]"], expression_terms(symbolics.add([])))

    def test_solver_only_edits_show_in_every_entry_that_holds_them(self):
        # None of these edits reaches cobra's own reaction data, so only the solver's side of the snapshot sees them.
        for solver in ("gurobi", "glpk"):
            with self.subTest(solver=solver):
                m = self.model(solver)
                r305, row = m.reactions.R305, m.constraints["m10[C_cy]"]
                forward = r305.forward_variable
                # The same bounds written in the other numeric type (0.0 for 0, 1000 for 1000.0) are no change.
                printed = str(forward)
                retyped = [float(bound) if isinstance(bound, int) else int(bound) for bound in (forward.lb, forward.ub)]
                self.assertEqual(self.changed(m, lambda: forward.set_bounds(*retyped)), (set(), set()))
                self.assertNotEqual(str(forward), printed)
                # A changed coefficient shows in its own row and nowhere else.
                self.assertEqual(self.changed(m, lambda: row.set_linear_coefficients({forward: 0.5})),
                                 ({"constraints"}, {row.name}))
                # A changed bound shows in every entry that holds the variable, also when that entry is compared alone.
                holding = {metabolite.id for metabolite in r305.metabolites} | {"one_sided", "two_sided"}
                self.assertEqual(self.changed(m, lambda: setattr(forward, "ub", 123)),
                                 ({"objective", "variables", "constraints"}, holding))
                # So does a changed type. (optlang then writes the bounds as ints, which alone is no change.)
                self.assertEqual(self.changed(m, lambda: setattr(forward, "type", "integer")),
                                 ({"objective", "variables", "constraints"}, holding))
                # With an integer variable in the model, the terms still say what the symbolic rows say.
                self.assertEqual(semantics(m)["objective"][1][0], (r305.id, 1.0, 0, 123, "integer"))
                self.assertEqual(constraint_terms(m), {c.name: expression_terms(c.expression) for c in m.constraints})


if __name__ == "__main__":
    unittest.main()
