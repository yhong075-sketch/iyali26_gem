"""Methylcitrate cycle curation on 7ccbb568: preservation, rejection, balance, SBML roundtrip, knockout logic, opt-in
flag; no LPs.

The baseline is the tracked candidate E5_353cefdc_xref_20261115.xml (7ccbb568). The checks are those of the other
line's test_methylcitrate_cycle.py, which builds its cbf2ef95 baseline through steps this branch does not have.
"""

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.methylcitrate_cycle import (CLOSED_ID, CURATION_PATH, GENES, GPR_SOURCE, NEW_ID,
                                                      NEW_STOICHIOMETRY, apply_methylcitrate_cycle)
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = "E5_353cefdc_xref_20261115.xml"
BASELINE_SHA256 = "7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146"
# The data file is copied byte-identical from the other line and records that line's baseline.
DATA_BASELINE_SHA256 = "cbf2ef95143bc5b74941257b9ddef38b5209f89f557cd591ac7ee75947934526"


def state(model):
    genes = {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}
    groups = {g.name: sorted(x.id for x in g.members) for g in model.groups}
    return semantics(model), annotations(model), genes, groups


def read_baseline():
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class MethylcitrateCycleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_baseline()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_baseline_file_and_the_inferred_rule(self):
        self.assertEqual(hashlib.sha256(MODEL.candidate_file(BASELINE).read_bytes()).hexdigest(), BASELINE_SHA256)
        self.assertEqual(self.spec["baseline_model_sha256"], DATA_BASELINE_SHA256)
        new = self.spec["new_reaction"]
        self.assertEqual(new["gpr_status"], "inferred")
        self.assertIn("INFERRED", new["notes"]["methylcitrate_cycle_gpr_status"])
        self.assertEqual(new["gpr"], self.model.reactions.get_by_id(GPR_SOURCE).gene_reaction_rule)
        self.assertIn("INFERRED", self.spec["closed_reaction"]["notes"]["methylcitrate_cycle_closure_evidence"])

    def test_applies_only_the_two_targets_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, before = read_baseline(), state(self.model)
        result = apply_methylcitrate_cycle(model)
        self.assertEqual(result["status"], "applied")
        after = state(model)
        added = model.reactions.get_by_id(NEW_ID)
        self.assertEqual({m.id: c for m, c in added.metabolites.items()}, NEW_STOICHIOMETRY)
        self.assertEqual(added.bounds, (0.0, 1000.0))
        self.assertEqual(added.gene_reaction_rule, model.reactions.get_by_id(GPR_SOURCE).gene_reaction_rule)
        self.assertEqual(model.reactions.get_by_id(CLOSED_ID).bounds, (0.0, 0.0))
        # Everything else is untouched: reactions other than the targets, all species, other genes and groups.
        sem_b, ann_b, genes_b, groups_b = before
        sem_a, ann_a, genes_a, groups_a = after
        self.assertEqual(sem_a["species"], sem_b["species"])
        self.assertEqual({k: v for k, v in sem_a["reactions"].items() if k not in (NEW_ID, CLOSED_ID)},
                         {k: v for k, v in sem_b["reactions"].items() if k != CLOSED_ID})
        self.assertEqual({k: v for k, v in ann_a.items() if k not in (NEW_ID, CLOSED_ID)},
                         {k: v for k, v in ann_b.items() if k != CLOSED_ID})
        self.assertEqual(ann_a[CLOSED_ID][1], ann_b[CLOSED_ID][1])  # annotation unchanged
        self.assertEqual({k: v for k, v in genes_a.items() if k not in GENES},
                         {k: v for k, v in genes_b.items() if k not in GENES})
        for gene in GENES:
            self.assertEqual(genes_a[gene][1], sorted(genes_b[gene][1] + [NEW_ID]))
        self.assertEqual({k: v for k, v in groups_a.items() if k != "Propanoate metabolism"},
                         {k: v for k, v in groups_b.items() if k != "Propanoate metabolism"})
        self.assertEqual(groups_a["Propanoate metabolism"], sorted(groups_b["Propanoate metabolism"] + [NEW_ID]))
        self.assertEqual([g.id for g in model.genes], [g.id for g in self.model.genes])
        self.assertEqual(apply_methylcitrate_cycle(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mcc.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == after, "Unexpected SBML roundtrip change")
        self.assertEqual(apply_methylcitrate_cycle(loaded)["status"], "already_correct")

    def test_new_reaction_is_mass_and_charge_balanced(self):
        model = self.model.copy()
        apply_methylcitrate_cycle(model)
        self.assertEqual(model.reactions.get_by_id(NEW_ID).check_mass_balance(), {})
        with model:
            model.reactions.get_by_id(NEW_ID).add_metabolites({model.metabolites.get_by_id("m28[C_mi]"): -1})
            self.assertEqual(model.reactions.get_by_id(NEW_ID).check_mass_balance(), {"charge": -1.0, "H": -1.0})

    def test_stale_inputs_fail_closed_without_edits(self):
        expected = {
            "formula": r"m78\[C_mi\] identity differs",
            "charge": r"m73\[C_mi\] identity differs",
            "missing_species": r"m42\[C_mi\] missing",
            "pdh_rule": "R693 is not the pyruvate dehydrogenase",
            "pdh_stoichiometry": "R693 is not the pyruvate dehydrogenase",
            "r2075_bounds": "R2075 bounds or notes differ",
            "r2075_equation": "R2075 methylcitrate cycle precondition differs",
            "r2075_gene": "R2075 methylcitrate cycle precondition differs",
            "r2075_closed_elsewhere": "R2075 is already closed without this curation",
            "r2075_extra_note": "R2075 bounds or notes differ",
            "context": r"m42\[C_mi\] reaction context differs",
            "context_cy": r"m2047\[C_cy\] reaction context differs",
            "duplicate": r"Reactions with the new chemistry already exist \['DUP'\]",
            "duplicate_reverse": r"Reactions with the new chemistry already exist \['DUP'\]",
            "new_other_definition": "R_OBDHm already exists with another definition",
            "partial_new": "Only part of the methylcitrate cycle curation is present",
            "partial_closed": "Only part of the methylcitrate cycle curation is present",
            "gene_context": "YALI1E31947g reaction context differs",
            "group": "Group 'Propanoate metabolism' is missing or not unique",
        }
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            m = bad.metabolites.get_by_id
            if field == "formula":
                m("m78[C_mi]").formula = "C24H40N7O17P3S"
            elif field == "charge":
                m("m73[C_mi]").charge = -4
            elif field == "missing_species":
                bad.remove_metabolites([m("m42[C_mi]")])
            elif field == "pdh_rule":
                bad.reactions.R693.gene_reaction_rule = "YALI1F27556g"
            elif field == "pdh_stoichiometry":
                bad.reactions.R693.add_metabolites({m("m28[C_mi]"): 1})
            elif field == "r2075_bounds":
                bad.reactions.R2075.bounds = (-1000.0, 1000.0)
            elif field == "r2075_equation":
                bad.reactions.R2075.add_metabolites({m("m28[C_mi]"): 1})
            elif field == "r2075_gene":
                bad.reactions.R2075.gene_reaction_rule = "YALI1F27556g"
            elif field == "r2075_closed_elsewhere":
                bad.reactions.R2075.bounds = (0.0, 0.0)
            elif field == "r2075_extra_note":
                bad.reactions.R2075.notes = {"unreviewed_note": "added elsewhere"}
            elif field == "context":
                extra = bad.reactions.R515.copy()
                extra.id = "EXTRA"
                bad.add_reactions([extra])
            elif field == "context_cy":
                extra = bad.reactions.R2075.copy()
                extra.id = "EXTRA"
                bad.add_reactions([extra])
            elif field in ("duplicate", "duplicate_reverse"):
                from cobra import Reaction
                dup = Reaction("DUP")
                bad.add_reactions([dup])
                sign = 1 if field == "duplicate" else -1
                dup.add_metabolites({m(k): sign * v for k, v in NEW_STOICHIOMETRY.items() if k != "m28[C_mi]"})
            elif field == "new_other_definition":
                apply_methylcitrate_cycle(bad)
                bad.reactions.get_by_id(NEW_ID).gene_reaction_rule = "YALI1F27556g"
            elif field == "partial_new":
                apply_methylcitrate_cycle(bad)
                bad.reactions.R2075.bounds = (0.0, 1000.0)
                bad.reactions.R2075.notes = {}
            elif field == "partial_closed":
                apply_methylcitrate_cycle(bad)
                bad.remove_reactions([bad.reactions.get_by_id(NEW_ID)])
            elif field == "gene_context":
                bad.reactions.R13.gene_reaction_rule = f"{bad.reactions.R13.gene_reaction_rule} or YALI1E31947g"
            else:
                group = next(g for g in bad.groups if g.name == "Propanoate metabolism")
                bad.remove_groups([group])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_methylcitrate_cycle(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        edits = {
            "status": lambda s: s.update(status="draft"),
            "schema_version": lambda s: s.update(schema_version=2),
            "curation_id": lambda s: s.update(curation_id="METHYLCITRATE"),
            "new_id": lambda s: s["new_reaction"].update(reaction_id="R_OTHER"),
            "bounds": lambda s: s["new_reaction"].update(bounds=[-1000.0, 1000.0]),
            "stoichiometry": lambda s: s["new_reaction"]["stoichiometry"].update({"m28[C_mi]": 3.0}),
            "gpr_status": lambda s: s["new_reaction"].update(gpr_status="verified"),
            "gpr_note": lambda s: s["new_reaction"]["notes"].update(methylcitrate_cycle_gpr_status="verified"),
            "gpr": lambda s: s["new_reaction"].update(gpr="YALI1F27556g"),
            "closed_after": lambda s: s["closed_reaction"].update(after_bounds=[0.0, 1.0]),
            "closed_before": lambda s: s["closed_reaction"].update(before_bounds=[-1000.0, 1000.0]),
            "previous_bounds_note": lambda s: s["closed_reaction"]["notes"].update(
                methylcitrate_cycle_previous_bounds="[0, 0]"),
            "note_angle": lambda s: s["new_reaction"]["notes"].update(methylcitrate_cycle_protons="a -> b"),
            "note_apostrophe": lambda s: s["closed_reaction"]["notes"].update(
                methylcitrate_cycle_curation="the model's route"),
            "empty_note": lambda s: s["new_reaction"]["notes"].update(methylcitrate_cycle_evidence=""),
            "genes": lambda s: s["genes"].pop(GENES[0]),
            "gene_after": lambda s: s["genes"][GENES[1]].update(reactions_after=["R693"]),
            "species": lambda s: s["species"].pop("m44[C_mi]"),
            "context": lambda s: s["species_reactions_before"].pop("m78[C_mi]"),
        }
        for key, edit in edits.items():
            spec = copy.deepcopy(self.spec)
            edit(spec)
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_methylcitrate_cycle(self.model.copy(), spec)

    def test_each_pdh_gene_knockout_closes_r693_and_the_new_step(self):
        model = self.model.copy()
        apply_methylcitrate_cycle(model)
        for gene_id in GENES:
            with model, self.subTest(gene=gene_id):
                model.genes.get_by_id(gene_id).knock_out()
                self.assertEqual(model.reactions.get_by_id(GPR_SOURCE).bounds, (0, 0))
                self.assertEqual(model.reactions.get_by_id(NEW_ID).bounds, (0, 0))

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_methylcitrate_cycle(model)["status"], "applied")
            self.assertIn(NEW_ID, model.reactions)
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertEqual(apply_methylcitrate_cycle(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).methylcitrate_cycle)
        self.assertTrue(parse_args(["--methylcitrate-cycle"]).methylcitrate_cycle)


if __name__ == "__main__":
    unittest.main()
