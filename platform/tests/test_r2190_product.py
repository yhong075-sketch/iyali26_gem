"""In-memory checks for the 2026-10-19 R2190 product correction (R2190-RIBULOSE); no optimization."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.coq9 import exact_residual
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.r2190_product import SPEC_PATH, apply_r2190_product
from scripts.gem_annotate.sbml import write_deterministic_sbml_model

SPEC = json.loads(SPEC_PATH.read_text())
NEW, OLD = SPEC["new_species"]["id"], SPEC["replaced_species"]["id"]


def _reaction(model, rid, stoichiometry, bounds=(0.0, 1000.0), gpr=""):
    reaction = Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([reaction])
    reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in stoichiometry.items()})
    reaction.gene_reaction_rule = gpr
    return reaction


def fixture():
    """R2190 in its a56f3649 state plus one stand-in for each other reaction of D-arabinose."""
    model = Model("r2190_fixture")
    model.add_metabolites([Metabolite(mid, formula=p["formula"], charge=p["charge"], compartment=p["compartment"])
                           for mid, p in SPEC["species"].items()])
    model.metabolites.get_by_id(OLD).name = SPEC["replaced_species"]["name"]
    reaction = _reaction(model, "R2190", SPEC["before_stoichiometry"], SPEC["bounds"], SPEC["gpr"])
    reaction.name = SPEC["reaction_name"]
    reaction.annotation = {"ec-code": SPEC["ec_code"], "kegg.reaction": "R07134"}
    for rid in SPEC["replaced_species"]["reactions_after"]:
        _reaction(model, rid, {OLD: -1.0})
    return model


def definition(model):
    return model_definition(model), {m.id: (m.annotation, m.notes) for m in model.metabolites}, \
        {r.id: (r.annotation, r.notes, r.gene_reaction_rule) for r in model.reactions}


class R2190ProductTests(unittest.TestCase):
    def test_replaces_the_product_only_and_is_idempotent(self):
        model = fixture()
        before = model_definition(model)
        result = apply_r2190_product(model)
        self.assertEqual(result["status"], "applied")
        after = model_definition(model)
        self.assertEqual({k: v for k, v in after["reactions"].items() if k != "R2190"},
                         {k: v for k, v in before["reactions"].items() if k != "R2190"})
        self.assertEqual({k: v for k, v in after["metabolites"].items() if k != NEW}, before["metabolites"])
        reaction = model.reactions.R2190
        self.assertEqual({m.id: c for m, c in reaction.metabolites.items()}, SPEC["after_stoichiometry"])
        self.assertEqual((reaction.bounds, reaction.gene_reaction_rule, reaction.name),
                         (tuple(SPEC["bounds"]), SPEC["gpr"], SPEC["reaction_name"]))
        self.assertEqual(exact_residual(reaction), {})
        self.assertEqual({k: reaction.notes[k] for k in SPEC["notes"]}, SPEC["notes"])
        met = model.metabolites.get_by_id(NEW)
        for field in ("name", "formula", "charge", "compartment", "annotation", "notes"):
            self.assertEqual(getattr(met, field), SPEC["new_species"][field])
        self.assertEqual(sorted(r.id for r in met.reactions), ["R2190"])
        self.assertEqual(result["replaced_species_usage_after"], SPEC["replaced_species"]["reactions_after"])
        state = definition(model)
        self.assertEqual(apply_r2190_product(model)["status"], "already_correct")
        self.assertEqual(definition(model), state)

    def test_sbml_roundtrip_keeps_the_after_state(self):
        model = fixture()
        apply_r2190_product(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r2190.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertEqual(model_definition(loaded), model_definition(model))
        self.assertEqual(loaded.metabolites.get_by_id(NEW).annotation, SPEC["new_species"]["annotation"])
        self.assertEqual(apply_r2190_product(loaded)["status"], "already_correct")

    def test_stale_inputs_fail_closed_without_edits(self):
        def name(m): m.reactions.R2190.name = "D-arabinitol dehydrogenase"
        def gpr(m): m.reactions.R2190.gene_reaction_rule = "YALI1E15452g"
        def ec(m): m.reactions.R2190.annotation["ec-code"] = "1.1.1.11"
        def bounds(m): m.reactions.R2190.bounds = (0.0, 1000.0)
        def proton(m):
            m.add_metabolites([Metabolite("m10[C_cy]", formula="H", charge=1, compartment="C_cy")])
            m.reactions.R2190.add_metabolites({m.metabolites.get_by_id("m10[C_cy]"): 1.0})
        def formula(m): m.metabolites.get_by_id(OLD).formula = "C5H11O8P"
        def replaced_name(m): m.metabolites.get_by_id(OLD).name = "D-ribulose_C5H10O5"
        def other_user(m): _reaction(m, "extra", {OLD: 1.0})
        def id_taken(m): m.add_metabolites([Metabolite("m2061[C_mi]", compartment="C_mi")])
        def note(m): m.reactions.R2190.notes["r2190_product_curation"] = "another review"
        definition_differs = "R2190 definition differs"
        cases = ((name, definition_differs), (gpr, definition_differs), (ec, definition_differs),
                 (bounds, definition_differs), (proton, "stoichiometry differs from the curated before-state"),
                 (formula, r"m454\[C_cy\]: species differs"), (replaced_name, r"m454\[C_cy\]: name differs"),
                 (other_user, r"m454\[C_cy\]: usage differs"), (id_taken, "species ID already in use"),
                 (note, "Conflicting R2190 product notes"))
        for change, message in cases:
            with self.subTest(change.__name__):
                model = fixture()
                change(model)
                before = definition(model)
                with self.assertRaisesRegex(ValueError, message):
                    apply_r2190_product(model)
                self.assertEqual(definition(model), before)

    def test_unexpected_or_unbalanced_curation_is_rejected(self):
        scope = copy.deepcopy(SPEC)
        scope["reaction_id"] = "R1964"
        extra_product = copy.deepcopy(SPEC)
        extra_product["after_stoichiometry"]["m10[C_cy]"] = 1.0
        changed_coefficient = copy.deepcopy(SPEC)
        changed_coefficient["after_stoichiometry"]["m123[C_cy]"] = 2.0
        unbalanced = copy.deepcopy(SPEC)
        unbalanced["new_species"]["formula"] = "C5H12O5"
        cases = ((scope, "Unexpected R2190 product curation scope"),
                 (extra_product, "Unexpected R2190 product curation scope"),
                 (changed_coefficient, "Unexpected R2190 product curation scope"),
                 (unbalanced, "not element and charge balanced"))
        for spec, message in cases:
            model = fixture()
            before = definition(model)
            with self.assertRaisesRegex(ValueError, message):
                apply_r2190_product(model, spec)
            self.assertEqual(definition(model), before)

    def test_rerun_rejects_a_drifted_after_state(self):
        model = fixture()
        apply_r2190_product(model)
        model.add_metabolites([Metabolite("m2061[C_mi]", compartment="C_mi")])
        with self.assertRaisesRegex(ValueError, "after-state differs"):
            apply_r2190_product(model)

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r2190_ribulose_product)
        self.assertTrue(parse_args(["--r2190-ribulose-product"]).r2190_ribulose_product)


if __name__ == "__main__":
    unittest.main()
