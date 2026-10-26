"""Thiamine phosphate chemistry: preservation, rejection, SBML roundtrip, exact edits, balance and the flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.thiamine_phosphate_chemistry import CURATION_PATH, apply_thiamine_phosphate_chemistry
from tests.test_coq9_curation import annotations, semantics
from tests.test_gpi_chain_links import group_table

BASELINE = "E5_coq9_resp_pkg_r132_combined_gpi_tpp_v2_20261022.xml"


def state(model):
    species = {m.id: (copy.deepcopy(m.annotation), copy.deepcopy(m.notes)) for m in model.metabolites}
    return semantics(model), annotations(model), group_table(model), species


def balance(model, rid):
    reaction = model.reactions.get_by_id(rid)
    return {k: v for k, v in reaction.check_mass_balance().items() if abs(v) > 1e-9}


class ThiaminePhosphateChemistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        precondition = "precondition differs"
        expected = {"no_tpp_requirement": "requires the TPP requirement",
                    "context_stoichiometry": "R752 differs from the reviewed reaction context",
                    "context_new_reaction": "m643\\[C_cy\\] identity or reaction context differs",
                    "species_name": "m977\\[C_cy\\] identity or reaction context differs",
                    "species_formula": f"m975\\[C_cy\\] {precondition}",
                    "species_annotation": f"m977\\[C_cy\\] {precondition}",
                    "r1342_annotation": f"R1342 {precondition}",
                    "extra_note": "m643\\[C_cy\\] notes outside this curation differ",
                    "conflicting_note": "Conflicting R1342 thiamine phosphate chemistry notes",
                    "after_without_notes": "m976\\[C_cy\\] already carries the change without this curation",
                    "notes_without_change": "m971\\[C_ex\\] carries this curation's notes without the change",
                    "partial": "Only part of the thiamine phosphate chemistry is present",
                    "partner_migrated": "m141\\[C_cy\\] on R\\d+ differs from the reviewed partner chemistry"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            if field == "no_tpp_requirement":
                tpp = bad.metabolites.get_by_id("m975[C_cy]")
                bad.reactions.biomass_C.subtract_metabolites({tpp: bad.reactions.biomass_C.get_coefficient(tpp)})
            elif field == "context_stoichiometry":
                # Removing R752's H+ was considered and not approved.
                bad.reactions.R752.subtract_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "context_new_reaction":
                source = Reaction("R_TEST_SOURCE", lower_bound=0.0, upper_bound=1000.0)
                bad.add_reactions([source])
                source.add_metabolites({bad.metabolites.get_by_id("m643[C_cy]"): 1})
            elif field == "species_name":
                bad.metabolites.get_by_id("m977[C_cy]").name = "thiamine monophosphate"
            elif field == "species_formula":
                # The pH 7.3 form, an option that was not chosen.
                tpp = bad.metabolites.get_by_id("m975[C_cy]")
                tpp.formula, tpp.charge = "C12H16N4O7P2S", -2
            elif field == "species_annotation":
                tmp = bad.metabolites.get_by_id("m977[C_cy]")
                tmp.annotation = {**tmp.annotation, "kegg.compound": "C00363"}
            elif field == "r1342_annotation":
                bad.reactions.R1342.annotation = {**bad.reactions.R1342.annotation, "ec-code": "2.7.4.9"}
            elif field == "extra_note":
                bad.metabolites.get_by_id("m643[C_cy]").notes = {"unreviewed_note": "added elsewhere"}
            elif field == "conflicting_note":
                bad.reactions.R1342.notes = {**bad.reactions.R1342.notes,
                                             "thiamine_phosphate_chemistry_status": "another review"}
            elif field == "after_without_notes":
                ttp, after = bad.metabolites.get_by_id("m976[C_cy]"), self.spec["species"]["m976[C_cy]"]["after"]
                ttp.formula, ttp.charge, ttp.annotation = after["formula"], after["charge"], dict(after["annotation"])
            elif field == "notes_without_change":
                apply_thiamine_phosphate_chemistry(bad)
                tdp, before = bad.metabolites.get_by_id("m971[C_ex]"), self.spec["species"]["m971[C_ex]"]["before"]
                tdp.formula, tdp.charge, tdp.annotation = before["formula"], before["charge"], dict(before["annotation"])
            elif field == "partner_migrated":
                # ATP at its pH 7.3 form: the KEGG forms would no longer balance, so the review is stale.
                atp = bad.metabolites.get_by_id("m141[C_cy]")
                atp.formula, atp.charge = "C10H12N5O13P3", -4
            else:
                apply_thiamine_phosphate_chemistry(bad)
                hetp, before = bad.metabolites.get_by_id("m643[C_cy]"), self.spec["species"]["m643[C_cy]"]["before"]
                hetp.formula, hetp.charge, hetp.annotation = before["formula"], before["charge"], dict(before["annotation"])
                hetp.notes = {}
            before_state = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_thiamine_phosphate_chemistry(bad)
            self.assertEqual(state(bad), before_state)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for key, value in (("status", "user_authorized"), ("curation_id", "THIAMINE-P-CHEMISTRY"),
                           ("species", {"m977[C_cy]": self.spec["species"]["m977[C_cy]"]}),
                           ("reactions", {}), ("context_reactions", {}), ("schema_version", 1.0),
                           ("partner_species", None)):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_thiamine_phosphate_chemistry(self.model.copy(), spec)

    def test_spec_values_outside_the_reviewed_change_are_rejected_without_edits(self):
        def charge_not_integer(spec):
            spec["species"]["m977[C_cy]"]["after"]["charge"] = 1.5

        def missing_kegg(spec):
            del spec["species"]["m643[C_cy]"]["after"]["annotation"]["kegg.compound"]

        def after_equals_before(spec):
            row = spec["species"]["m1116[C_mi]"]
            row["after"] = copy.deepcopy(row["before"])

        def compartments_disagree(spec):
            # The same molecule must get the same chemistry in every compartment.
            spec["species"]["m971[C_ex]"]["after"]["charge"] = -2

        def r1342_other_annotation(spec):
            spec["reactions"]["R1342"]["annotation_after"]["kegg.reaction"] = "R00616"

        def dtmp_formula(spec):
            # A consistent but wrong identity in both compartments.
            for mid in ("m977[C_cy]", "m973[C_ex]"):
                after = spec["species"][mid]["after"]
                after["formula"], after["charge"] = "C10H15N2O8P", 0
                after["annotation"]["kegg.compound"] = "C00364"

        def hetp_ph73(spec):
            after = spec["species"]["m643[C_cy]"]["after"]
            after["formula"], after["charge"] = "C6H8NO4PS", -2

        def molecule_relabelled(spec):
            spec["species"]["m971[C_ex]"]["molecule"] = "thiamine diphosphate, extracellular"
            spec["species"]["m971[C_ex]"]["after"]["charge"] = -2

        def charge_true(spec):
            spec["species"]["m643[C_cy]"]["after"]["charge"] = True

        model = self.model.copy()
        before = state(model)
        for edit in (charge_not_integer, missing_kegg, after_equals_before, compartments_disagree,
                     r1342_other_annotation, dtmp_formula, hetp_ph73, molecule_relabelled, charge_true):
            spec = copy.deepcopy(self.spec)
            edit(spec)
            with self.subTest(edit=edit.__name__), self.assertRaisesRegex(ValueError, "spec values differ"):
                apply_thiamine_phosphate_chemistry(model, spec)
            self.assertEqual(state(model), before)

    def test_changes_only_the_targets_balances_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        for mid, row in self.spec["species"].items():
            metabolite = expected.metabolites.get_by_id(mid)
            metabolite.formula, metabolite.charge = row["after"]["formula"], row["after"]["charge"]
            metabolite.annotation = dict(row["after"]["annotation"])
            metabolite.notes = self._notes(row)
        r1342 = expected.reactions.R1342
        r1342.annotation = dict(self.spec["reactions"]["R1342"]["annotation_after"])
        r1342.notes = {**r1342.notes, **self._notes(self.spec["reactions"]["R1342"])}
        # Exchanges (R1339, R1341) and biomass_C are not mass-balanced by design.
        before_balance = {rid: balance(model, rid) for rid in self.spec["context_reactions"]
                          if rid != "biomass_C" and not model.reactions.get_by_id(rid).boundary}
        result = apply_thiamine_phosphate_chemistry(model)
        self.assertEqual((result["status"], result["changed_reactions"]), ("applied", ["R1342"]))
        self.assertEqual(result["changed_species"], sorted(self.spec["species"]))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        after_balance = {rid: balance(model, rid) for rid in before_balance}
        balanced = sorted(rid for rid, residual in after_balance.items() if not residual)
        self.assertEqual(balanced, sorted(self.spec["expected_balanced_after"]))
        # R752, R408 and R750 keep the one H+ of their pH 7.3 equations.
        for rid in ("R752", "R408", "R750"):
            self.assertEqual(after_balance[rid], {"charge": 1.0, "H": 1.0})
        for rid in before_balance:
            if not before_balance[rid]:
                self.assertFalse(after_balance[rid], f"{rid} became unbalanced")
        self.assertEqual(model.reactions.R1342.annotation["ec-code"], "2.7.4.16")
        self.assertEqual(apply_thiamine_phosphate_chemistry(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "thiamine_phosphate_chemistry.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_thiamine_phosphate_chemistry(loaded)["status"], "already_correct")

    def test_flag_requires_the_tpp_requirement(self):
        self.assertFalse(parse_args([]).thiamine_phosphate_chemistry)
        self.assertTrue(parse_args(["--thiamine-phosphate-chemistry"]).thiamine_phosphate_chemistry)
        with self.assertRaisesRegex(ValueError, "defined on the TPP requirement candidate"):
            build_reference_chain(thiamine_phosphate_chemistry=True)

    def _notes(self, row):
        return {"thiamine_phosphate_chemistry_curation_id": self.spec["curation_id"],
                "thiamine_phosphate_chemistry_status": self.spec["status"],
                "thiamine_phosphate_chemistry_authorization": self.spec["authorization"],
                "thiamine_phosphate_chemistry_evidence": row["evidence"],
                "thiamine_phosphate_chemistry_evidence_tier": row["evidence_tier"],
                "thiamine_phosphate_chemistry_rationale": row["rationale"],
                "thiamine_phosphate_chemistry_previous_state": row["previous_state"]}


if __name__ == "__main__":
    unittest.main()
