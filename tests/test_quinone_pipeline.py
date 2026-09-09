"""Integration checks for the existing Q9 chain and its annotation prerequisites."""

import copy
import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.core.gene import GPR

from scripts.gem_annotate.coq9 import boolean_key
from scripts.gem_annotate.metabolites import _apply_mnxm, fix_proton_water_balance
from scripts.gem_annotate.quinone import (
    _COQ9_ROUTE_IDS, apply_reviewed_quinone_step_gprs,
    correct_external_ndh2_gpr_and_remove_duplicate, remove_spurious_quinone_branches,
    replace_coq6_route_with_coq9, run_quinone_step,
)
from tests.test_coq9_curation import annotations, semantics

FIXTURE = Path(__file__).parent / "fixtures" / "quinone_legacy.json"


def legacy_model():
    data = json.loads(FIXTURE.read_text())
    model = Model("quinone_annotated_precursor")
    for mid, row in data["metabolites"].items():
        met = Metabolite(mid, name=row["name"], formula=row["formula"],
                         charge=row["charge"], compartment=row["compartment"])
        met.annotation = row["annotation"]
        model.add_metabolites([met])
    for rid, row in data["reactions"].items():
        reaction = Reaction(rid, name=row["name"])
        reaction.bounds = row["bounds"]
        reaction.add_metabolites({model.metabolites.get_by_id(mid): value
                                 for mid, value in row["stoichiometry"].items()})
        reaction.gene_reaction_rule = row["gpr"]
        reaction.annotation, reaction.notes = row["annotation"], row["notes"]
        model.add_reactions([reaction])
    model.objective = "R305"
    return model, data


class QuinonePipelineTests(unittest.TestCase):
    def test_legacy_balancer_compartment_ties_are_reproducible_across_processes(self):
        script = '''
import json
from cobra import Model, Metabolite, Reaction
from scripts.gem_annotate.metabolites import fix_proton_water_balance
m = Model('tie')
a = Metabolite('a', formula='CH', charge=0, compartment='C_mi')
b = Metabolite('b', formula='C', charge=0, compartment='C_cy')
for comp in ('C_mi', 'C_cy'):
 h = Metabolite('h_'+comp, formula='H', charge=1, compartment=comp)
 h.annotation['bigg.metabolite'] = 'h'
 m.add_metabolites([h])
r = Reaction('transport'); r.add_metabolites({a:-1,b:1}); m.add_reactions([r])
fix_proton_water_balance(m)
print(json.dumps({s.id:v for s,v in r.metabolites.items()},sort_keys=True))
'''
        results = []
        for seed in ('1', '2'):
            output = subprocess.check_output([sys.executable, '-c', script],
                                             env={**os.environ, 'PYTHONHASHSEED': seed}, text=True)
            results.append(json.loads(output.splitlines()[-1]))
        self.assertEqual(results[0], results[1])
        self.assertEqual(results[0]['h_C_cy'], 1)

    def test_existing_chain_reproduces_saved_local_mathematics_and_is_idempotent(self):
        model, data = legacy_model()
        genes = {g.id for g in model.genes}
        proton_before = dict(semantics(model)["reactions"]["R305"][0])
        fix_proton_water_balance(model)
        self.assertEqual(proton_before, semantics(model)["reactions"]["R305"][0])
        operations = (replace_coq6_route_with_coq9,
                      correct_external_ndh2_gpr_and_remove_duplicate,
                      remove_spurious_quinone_branches, apply_reviewed_quinone_step_gprs)
        for operation in operations:
            self.assertNotEqual(run_quinone_step(model, operation)["status"], "conflict")
        self.assertEqual(genes, {g.id for g in model.genes})
        for rid, expected in data["expected"].items():
            reaction = model.reactions.get_by_id(rid)
            self.assertEqual({m.id: v for m, v in reaction.metabolites.items()}, expected["stoichiometry"], rid)
            self.assertEqual(list(reaction.bounds), expected["bounds"], rid)
            self.assertEqual(boolean_key(reaction.gpr.body), boolean_key(GPR.from_string(expected["gpr"]).body), rid)
            if rid in _COQ9_ROUTE_IDS:
                self.assertEqual(reaction.check_mass_balance(), {}, rid)
        state = (semantics(model), annotations(model))
        for operation in operations:
            self.assertEqual(run_quinone_step(model, operation)["status"], "already_correct")
        self.assertEqual(state, (semantics(model), annotations(model)))

    def test_changed_duplicate_and_branch_leave_input_untouched(self):
        for rid, operation in (("R2063", correct_external_ndh2_gpr_and_remove_duplicate),
                               ("R189", remove_spurious_quinone_branches)):
            model, _ = legacy_model()
            model.reactions.get_by_id(rid).upper_bound = 17
            before = (semantics(model), annotations(model),
                      {g.id: copy.deepcopy(g.annotation) for g in model.genes})
            self.assertEqual(run_quinone_step(model, operation)["status"], "conflict")
            self.assertEqual(before, (semantics(model), annotations(model),
                                     {g.id: g.annotation for g in model.genes}))

    def test_unexpected_gpr_is_not_replaced_even_if_legacy_markers_remain(self):
        model, _ = legacy_model()
        reaction = model.reactions.R570
        reaction.gene_reaction_rule = "(" + reaction.gene_reaction_rule + ") or engineered_gene"
        before = annotations(model)
        self.assertEqual(run_quinone_step(model, correct_external_ndh2_gpr_and_remove_duplicate)["status"], "conflict")
        self.assertEqual(before, annotations(model))

    def test_atomic_q9_chemistry_preserves_source_nad_and_other_curation_scope(self):
        model, _ = legacy_model()
        nad = model.metabolites.get_by_id("m27[C_mi]")
        before = (nad.formula, nad.charge)
        properties = {"MNXM_test": {"formula": "C21H26N7O14P2", "charge": "-1",
                                   "inchi": "", "inchikey": ""}}
        _apply_mnxm(nad, "MNXM_test", {}, properties)
        self.assertEqual((nad.formula, nad.charge), before)
        self.assertIn("metanetx_formula_charge_conflict", nad.notes)
        # The existing ER contract expects its old annotation path; Q9 does
        # not silently migrate that separate source tuple.
        partner = Metabolite("m1439[C_em]", formula="C21H29N7O17P3", charge=0)
        properties["MNXM_test"].update(formula="C21H25N7O17P3", charge="-3")
        _apply_mnxm(partner, "MNXM_test", {}, properties)
        self.assertEqual((partner.formula, partner.charge), ("C21H29N7O17P3", -3))


if __name__ == "__main__":
    unittest.main()
