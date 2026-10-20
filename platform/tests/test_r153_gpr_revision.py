"""In-memory checks for the 2026-10-19 R153 GPR revision; no optimization."""

import copy
import json
import unittest

from cobra import Metabolite, Model, Reaction

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.r153_gpr_revision import SPEC_PATH, apply_r153_gpr_revision

SPEC = json.loads(SPEC_PATH.read_text())


def fixture():
    """Merged R153 as in the baseline, plus R152 sharing argininosuccinate with its own gene."""
    before = SPEC["before"]
    model = Model("r153_fixture")
    model.add_metabolites([Metabolite(mid, formula=row["formula"], charge=row["charge"],
                                      compartment=row["compartment"])
                           for mid, row in before["species"].items()])
    r153 = Reaction("R153", lower_bound=before["bounds"][0], upper_bound=before["bounds"][1])
    r152 = Reaction("R152", lower_bound=-1000.0, upper_bound=1000.0)
    model.add_reactions([r153, r152])
    r153.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in before["stoichiometry"].items()})
    r153.gene_reaction_rule = before["gpr"]
    r153.notes = {"PROTEIN_CLASS": "6.3.4.5", **before["required_notes"],
                  **{k: row["before"] for k, row in SPEC["notes_replaced"].items()}}
    r152.add_metabolites({model.metabolites.get_by_id("m312[C_cy]"): -1.0})
    r152.gene_reaction_rule = "YALI1D34865g"
    model.genes.get_by_id(SPEC["gene"]["id"]).annotation.update(SPEC["gene"]["annotation"])
    return model


def state(model):
    return {"reactions": {r.id: (r.gene_reaction_rule, dict(r.notes), tuple(r.bounds),
                                 {m.id: c for m, c in r.metabolites.items()}) for r in model.reactions},
            "metabolites": {m.id: (m.formula, m.charge, m.compartment) for m in model.metabolites},
            "genes": {g.id: (dict(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}}


def _set_charge(model):
    model.metabolites.get_by_id("m312[C_cy]").charge = -1


def _add_r2176(model):
    model.add_reactions([Reaction("R2176")])


class R153GprRevisionTests(unittest.TestCase):
    def test_clears_gpr_keeps_gene_and_is_idempotent(self):
        model = fixture()
        r152 = state(model)["reactions"]["R152"]
        result = apply_r153_gpr_revision(model)
        self.assertEqual(result["status"], "applied")
        r153 = model.reactions.R153
        self.assertEqual(r153.gene_reaction_rule, "")
        self.assertIn("YALI1D17462g", model.genes)  # The gene is kept, unlinked.
        self.assertEqual(len(model.genes.get_by_id("YALI1D17462g").reactions), 0)
        for key, row in SPEC["notes_replaced"].items():
            self.assertEqual(r153.notes[key], row["after"])
        for key, value in SPEC["notes_added"].items():
            self.assertEqual(r153.notes[key], value)
        self.assertEqual(r153.notes["PROTEIN_CLASS"], "6.3.4.5")
        self.assertEqual(r153.notes["reaction_merge_id"], SPEC["before"]["required_notes"]["reaction_merge_id"])
        self.assertEqual(state(model)["reactions"]["R152"], r152)
        after = state(model)
        self.assertEqual(apply_r153_gpr_revision(model)["status"], "already_correct")
        self.assertEqual(state(model), after)

    def test_guards_fail_closed_without_edits(self):
        cases = {
            "gpr": lambda m: setattr(m.reactions.R153, "gene_reaction_rule", "YALIUNK2 or YALI1D17462g"),
            "bounds": lambda m: setattr(m.reactions.R153, "upper_bound", 999.0),
            "stoichiometry": lambda m: m.reactions.R153.add_metabolites({m.metabolites.get_by_id("m141[C_cy]"): -1.0}),
            "charge": _set_charge,
            "r2176_present": _add_r2176,
            "gene_refseq": lambda m: m.genes.get_by_id("YALI1D17462g").annotation.update(refseq="XP_000000.1"),
            "gene_on_other_reaction": lambda m: setattr(m.reactions.R152, "gene_reaction_rule",
                                                        "YALI1D34865g or YALI1D17462g"),
            "old_note_differs": lambda m: m.reactions.R153.notes.update(gpr_evidence_status="changed"),
            "merge_note_differs": lambda m: m.reactions.R153.notes.update(reaction_merge_id="other"),
            "new_note_present": lambda m: m.reactions.R153.notes.update(r153_gpr_revision="other"),
        }
        for name, mutate in cases.items():
            model = fixture()
            mutate(model)
            before = state(model)
            with self.subTest(case=name), self.assertRaises(ValueError):
                apply_r153_gpr_revision(model)
            self.assertEqual(state(model), before)
        model, spec = fixture(), copy.deepcopy(SPEC)
        spec["change_id"] = "OTHER"
        before = state(model)
        with self.subTest(case="spec_scope"), self.assertRaises(ValueError):
            apply_r153_gpr_revision(model, spec)
        self.assertEqual(state(model), before)

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).remove_r153_peptidase_gpr)
        self.assertTrue(parse_args(["--remove-r153-peptidase-gpr"]).remove_r153_peptidase_gpr)


if __name__ == "__main__":
    unittest.main()
