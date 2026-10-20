"""In-memory checks for the 2026-10-19 restore of four selection-cleared GPRs; no optimization."""

import copy
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
    SPEC_PATH, TARGETS, apply_selection_gpr_restore, restore_after_selection)

SPEC = json.loads(SPEC_PATH.read_text())
BASELINE = "E5_coq9_resp_pkg_r132_20261012.xml"


def fixture(state="cleared"):
    model = Model("selection_gpr_restore_fixture")
    model.add_metabolites([Metabolite(mid, formula=row["formula"], charge=row["charge"],
                                      compartment=row["compartment"])
                           for mid, row in SPEC["species"].items()])
    for rid in TARGETS:
        target = SPEC["reactions"][rid][state]
        reaction = Reaction(rid, lower_bound=target["bounds"][0], upper_bound=target["bounds"][1])
        model.add_reactions([reaction])
        reaction.add_metabolites({model.metabolites.get_by_id(mid): c
                                  for mid, c in target["stoichiometry"].items()})
        reaction.gene_reaction_rule = target["gpr"]
    for rid in TARGETS:
        gene = SPEC["reactions"][rid]["gene"]
        if gene["id"] not in model.genes:
            # Cleared genes have no reaction; keep them in the model as the baseline does.
            holder = Reaction(f"holder_{gene['id']}")
            model.add_reactions([holder])
            holder.gene_reaction_rule = gene["id"]
        g = model.genes.get_by_id(gene["id"])
        g.name = gene.get("name", g.name)
        g.annotation.update({k: v for k, v in gene.items() if k not in ("id", "name")})
    return model


def notes(model):
    return {rid: dict(model.reactions.get_by_id(rid).notes) for rid in TARGETS}


class SelectionGprRestoreTests(unittest.TestCase):
    def test_restores_balanced_targets_and_is_idempotent(self):
        model = fixture()
        genes, metabolites = {g.id for g in model.genes}, {m.id for m in model.metabolites}
        result = apply_selection_gpr_restore(model)
        self.assertEqual(result["status"], "applied")
        self.assertEqual([r["status"] for r in result["reactions"]], ["applied"] * 4)
        for rid in TARGETS:
            reaction, target = model.reactions.get_by_id(rid), SPEC["reactions"][rid]["restored"]
            self.assertEqual({m.id: c for m, c in reaction.metabolites.items()}, target["stoichiometry"])
            self.assertEqual(reaction.gene_reaction_rule, target["gpr"])
            self.assertEqual(reaction.check_mass_balance(), {})
            self.assertEqual(reaction.notes["selection_gpr_restore"], SPEC["change_id"])
        self.assertNotIn(model.metabolites.get_by_id("m10[C_cy]"), model.reactions.R1303.metabolites)
        self.assertEqual(model.reactions.R2274.metabolites[model.metabolites.get_by_id("m214[C_ex]")], -3.0)
        self.assertEqual({g.id for g in model.genes}, genes)  # No gene added.
        self.assertEqual({m.id for m in model.metabolites}, metabolites)  # No metabolite added.
        with model:
            model.genes.YALI1A08512g.knock_out()
            self.assertEqual(model.reactions.R1302.bounds, (0, 0))
            self.assertEqual(model.reactions.R1303.bounds, (0, 0))
        definition, state = model_definition(model), notes(model)
        rerun = apply_selection_gpr_restore(model)
        self.assertEqual(rerun["status"], "already_correct")
        self.assertEqual((model_definition(model), notes(model)), (definition, state))

    def test_rerun_after_sbml_round_trip_is_a_no_op(self):
        model = fixture()
        apply_selection_gpr_restore(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "restored.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        definition, state = model_definition(loaded), notes(loaded)
        self.assertEqual(apply_selection_gpr_restore(loaded)["status"], "already_correct")
        self.assertEqual((model_definition(loaded), notes(loaded)), (definition, state))

    def test_stale_inputs_fail_closed_without_edits(self):
        def bounds(m): m.reactions.R612.upper_bound = 999.0
        def other_gpr(m): m.reactions.R1302.gene_reaction_rule = "YALI1E31685g"
        def mixed_state(m): m.reactions.R2274.gene_reaction_rule = "YALI1A21372g"
        def extra_species(m): m.reactions.R1303.add_metabolites({m.metabolites.get_by_id("m82[C_cy]"): 1.0})
        def formula(m): m.metabolites.get_by_id("m920[C_cy]").formula = "C8H12N2O5P"
        def charge(m): m.metabolites.get_by_id("m852[C_cy]").charge = -2
        def gene_name(m): m.genes.YALI1A21372g.name = "LIP7"
        def gene_identity(m): m.genes.YALI1E31685g.annotation["uniprot"] = "conflicting_accession"
        def note(m): m.reactions.R612.notes["selection_gpr_restore"] = "other change"
        def missing(m): m.remove_reactions([m.reactions.R1302])
        for change in (bounds, other_gpr, mixed_state, extra_species, formula, charge,
                       gene_name, gene_identity, note, missing):
            with self.subTest(change.__name__):
                model = fixture()
                change(model)
                before = (model_definition(model),
                          {r.id: dict(r.notes) for r in model.reactions})
                with self.assertRaises(ValueError):
                    apply_selection_gpr_restore(model)
                self.assertEqual((model_definition(model),
                                  {r.id: dict(r.notes) for r in model.reactions}), before)

    def test_unexpected_scope_or_unbalanced_target_is_rejected(self):
        bad_id = copy.deepcopy(SPEC); bad_id["change_id"] = "OTHER"
        missing = copy.deepcopy(SPEC); del missing["reactions"]["R1302"]
        unbalanced = copy.deepcopy(SPEC)
        unbalanced["reactions"]["R2274"]["restored"]["stoichiometry"]["m214[C_ex]"] = -2.0
        for name, spec in (("change_id", bad_id), ("scope", missing), ("unbalanced", unbalanced)):
            with self.subTest(name):
                model = fixture()
                before = model_definition(model)
                with self.assertRaises(ValueError):
                    apply_selection_gpr_restore(model, spec)
                self.assertEqual(model_definition(model), before)

    def test_flag_is_opt_in_and_needs_a_complete_selection(self):
        self.assertFalse(parse_args([]).restore_selection_cleared_gprs)
        self.assertTrue(parse_args(["--restore-selection-cleared-gprs"]).restore_selection_cleared_gprs)
        model = fixture()
        before = model_definition(model)
        with self.assertRaises(ValueError):
            restore_after_selection(model, {"complete": False})
        self.assertEqual(model_definition(model), before)
        self.assertEqual(restore_after_selection(model, {"complete": True})["status"], "applied")

    def test_baseline_candidate_changes_only_the_targets(self):
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        before = model_definition(model)
        self.assertEqual(apply_selection_gpr_restore(model)["status"], "applied")
        after = model_definition(model)
        self.assertEqual(after["metabolites"], before["metabolites"])
        self.assertEqual({k: v for k, v in after.items() if k != "reactions"},
                         {k: v for k, v in before.items() if k != "reactions"})
        changed = sorted(r for r in after["reactions"] if after["reactions"][r] != before["reactions"][r])
        self.assertEqual(changed, sorted(TARGETS))
        for rid in TARGETS:
            self.assertEqual(model.reactions.get_by_id(rid).check_mass_balance(), {})


if __name__ == "__main__":
    unittest.main()
