"""In-memory checks for the R765/R766 change set that drops the transketolase N-terminal fragment
YALI1D02625g (2026-10-19); no optimization."""

import json
import tempfile
import unittest
from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.selection_gpr_restore import (
    CHANGE_SETS, TKL_FRAGMENT_GPR_RESTORE, apply_selection_gpr_restore, restore_after_selection)

BASELINE = "E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml"
TARGETS = ("R765", "R766")
SPEC = json.loads(CHANGE_SETS[TKL_FRAGMENT_GPR_RESTORE][0].read_text())


def fixture(state="cleared"):
    model = Model("tkl_fragment_fixture")
    model.add_metabolites([Metabolite(mid, formula=row["formula"], charge=row["charge"],
                                      compartment=row["compartment"])
                           for mid, row in SPEC["species"].items()])
    for rid, row in SPEC["reactions"].items():
        target = row[state]
        reaction = Reaction(rid, lower_bound=target["bounds"][0], upper_bound=target["bounds"][1])
        model.add_reactions([reaction])
        reaction.add_metabolites({model.metabolites.get_by_id(mid): c
                                  for mid, c in target["stoichiometry"].items()})
        reaction.gene_reaction_rule = target["gpr"]
    for row in SPEC["reactions"].values():
        for gene in row["genes"]:
            g = model.genes.get_by_id(gene["id"])
            g.name = gene["name"]
            g.annotation.update({k: v for k, v in gene.items() if k not in ("id", "name")})
    return model


def genes(model):
    return {g.id: (g.name, dict(g.annotation)) for g in model.genes}


def snapshot(model):
    return model_definition(model), {r.id: dict(r.notes) for r in model.reactions}, genes(model)


class TklFragmentRestoreTests(unittest.TestCase):
    def test_record_scope(self):
        self.assertEqual(SPEC["change_id"], TKL_FRAGMENT_GPR_RESTORE)
        self.assertEqual(tuple(sorted(SPEC["reactions"])), TARGETS)
        for rid, row in SPEC["reactions"].items():
            with self.subTest(rid):
                self.assertEqual(row["cleared"]["gpr"], "YALI1E07744g or YALI1D02625g")
                self.assertEqual(row["restored"]["gpr"], "YALI1E07744g")
                self.assertEqual(row["restored"]["stoichiometry"], row["cleared"]["stoichiometry"])
                self.assertEqual(row["restored"]["bounds"], row["cleared"]["bounds"])
                self.assertEqual(row["removed_genes"], ["YALI1D02625g"])
                self.assertEqual(row["evidence_tier"], "inferred")
                self.assertEqual(row["unchanged_equation_residual"], {"charge": 2.0})

    def test_removes_fragment_keeps_gene_and_is_idempotent(self):
        model = fixture()
        before_genes, metabolites = genes(model), {m.id for m in model.metabolites}
        result = apply_selection_gpr_restore(model, change_id=TKL_FRAGMENT_GPR_RESTORE)
        self.assertEqual(result["status"], "applied")
        for rid, row in SPEC["reactions"].items():
            reaction = model.reactions.get_by_id(rid)
            self.assertEqual(reaction.gene_reaction_rule, "YALI1E07744g")
            self.assertEqual({m.id: c for m, c in reaction.metabolites.items()}, row["cleared"]["stoichiometry"])
            self.assertEqual(reaction.check_mass_balance(), {"charge": 2.0})  # pre-existing, unchanged
            self.assertEqual(reaction.notes["selection_gpr_restore"], TKL_FRAGMENT_GPR_RESTORE)
        self.assertEqual(genes(model), before_genes)  # YALI1D02625g kept, no gene added
        self.assertEqual(model.genes.YALI1D02625g.reactions, frozenset())
        self.assertEqual({m.id for m in model.metabolites}, metabolites)
        state = snapshot(model)
        self.assertEqual(apply_selection_gpr_restore(model, change_id=TKL_FRAGMENT_GPR_RESTORE)["status"],
                         "already_correct")
        self.assertEqual(snapshot(model), state)

    def test_gpr_propagation(self):
        model = fixture()
        apply_selection_gpr_restore(model, change_id=TKL_FRAGMENT_GPR_RESTORE)
        for gene, closed in (("YALI1E07744g", True), ("YALI1D02625g", False)):
            with self.subTest(gene), model:
                model.genes.get_by_id(gene).knock_out()
                for rid in TARGETS:
                    self.assertEqual(model.reactions.get_by_id(rid).bounds == (0, 0), closed, rid)

    def test_rerun_after_sbml_round_trip_is_a_no_op(self):
        model = fixture()
        apply_selection_gpr_restore(model, change_id=TKL_FRAGMENT_GPR_RESTORE)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "restored.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        state = snapshot(loaded)
        self.assertEqual(apply_selection_gpr_restore(loaded, change_id=TKL_FRAGMENT_GPR_RESTORE)["status"],
                         "already_correct")
        self.assertEqual(snapshot(loaded), state)

    def test_stale_inputs_fail_closed_without_edits(self):
        def fragment_only(m): m.reactions.R765.gene_reaction_rule = "YALI1D02625g"
        def complex_rule(m): m.reactions.R766.gene_reaction_rule = "YALI1E07744g and YALI1D02625g"
        def equation(m): m.reactions.R765.add_metabolites({m.metabolites.get_by_id("m10[C_cy]"): -2.0})
        def bounds(m): m.reactions.R766.bounds = (0.0, 1000.0)
        def species_charge(m): m.metabolites.get_by_id("m938[C_cy]").charge = -2
        def gene_name(m): m.genes.YALI1E07744g.name = "TKL1"
        def note(m): m.reactions.R766.notes["selection_gpr_restore"] = "other change"
        def missing(m): m.remove_reactions([m.reactions.R766])
        neither = "neither the cleared nor the restored state"
        for change, message in ((fragment_only, neither), (complex_rule, neither), (equation, neither),
                                (bounds, neither), (species_charge, "species differs"),
                                (gene_name, "name differs"), (note, "conflicting restore notes"),
                                (missing, "R766: reaction missing")):
            with self.subTest(change.__name__):
                model = fixture()
                change(model)
                before = snapshot(model)
                with self.assertRaisesRegex(ValueError, message):
                    apply_selection_gpr_restore(model, change_id=TKL_FRAGMENT_GPR_RESTORE)
                self.assertEqual(snapshot(model), before)

    def test_flag_is_opt_in_and_needs_a_complete_selection(self):
        self.assertFalse(parse_args([]).remove_selection_tkl_fragment_gpr)
        self.assertTrue(parse_args(["--remove-selection-tkl-fragment-gpr"]).remove_selection_tkl_fragment_gpr)
        model = fixture()
        before = snapshot(model)
        with self.assertRaises(ValueError):
            restore_after_selection(model, {"complete": False}, TKL_FRAGMENT_GPR_RESTORE)
        self.assertEqual(snapshot(model), before)
        self.assertEqual(restore_after_selection(model, {"complete": True}, TKL_FRAGMENT_GPR_RESTORE)["status"],
                         "applied")

    def test_baseline_candidate_changes_only_the_targets(self):
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        before, before_genes = model_definition(model), genes(model)
        self.assertEqual(apply_selection_gpr_restore(model, change_id=TKL_FRAGMENT_GPR_RESTORE)["status"], "applied")
        after = model_definition(model)
        self.assertEqual({k: v for k, v in after.items() if k != "reactions"},
                         {k: v for k, v in before.items() if k != "reactions"})
        changed = sorted(r for r in after["reactions"] if after["reactions"][r] != before["reactions"][r])
        self.assertEqual(changed, list(TARGETS))
        self.assertEqual(genes(model), before_genes)
        self.assertEqual(model.genes.YALI1D02625g.reactions, frozenset())
        self.assertEqual(sorted(r.id for r in model.genes.YALI1E07744g.reactions), list(TARGETS))


if __name__ == "__main__":
    unittest.main()
