"""Bounded version selection, idempotency, round-trip and local-conflict checks."""

import html
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from cobra import Model, Metabolite, Reaction
from cobra.io import read_sbml_model, write_sbml_model

from scripts.gem_annotate.reaction_selection import (
    SELECTION_PATH, apply_metadata_reaction_selection, reaction_fields, _same,
)
from tests.test_coq9_curation import annotations, semantics


class ReactionSelectionTests(unittest.TestCase):
    def setUp(self):
        self.spec = json.loads(SELECTION_PATH.read_text())
        self.model = Model("reaction_selection")
        species = {mid: props for row in self.spec["reactions"].values() for mid, props in row["species"].items()}
        self.model.add_metabolites([Metabolite(mid, **props) for mid, props in species.items()])
        for rid, row in self.spec["reactions"].items():
            r = Reaction(rid)
            r.add_metabolites({self.model.metabolites.get_by_id(mid): c for mid, c in row["before"]["stoichiometry"].items()})
            r.bounds = row["before"]["bounds"]
            # Include target-only genes without changing the tested reaction membership.
            r.gene_reaction_rule = row["after"]["gpr"]
            self.model.add_reactions([r])
            r.gene_reaction_rule = row["before"]["gpr"]
            r.notes = {k: html.unescape(str(v)) for k, v in row["superseded_notes"].items()}
        self.model.reactions.biomass_C.notes = {"canonical_trna_biomass_mode": "split_v1",
                                              "experimental_trna_biomass_mode": "split_v1"}

    def test_exact_selection_preserves_other_fields_and_roundtrips(self):
        m = self.model
        m.objective = "biomass_C"
        m.add_cons_vars(m.problem.Constraint(m.reactions.R159.flux_expression, ub=17, name="extra_cap"))
        before = semantics(m)
        report = apply_metadata_reaction_selection(m)
        self.assertTrue(report["complete"])
        counts = Counter(f["field"] for row in report["items"] for f in row["fields"] if f["status"] == "applied")
        self.assertEqual(dict(counts), {"stoichiometry": 202, "bounds": 7, "gpr": 9})
        after = semantics(m)
        for field in ("species", "compartments", "genes", "objective"):
            self.assertEqual(before[field], after[field], field)
        self.assertEqual(before["constraints"]["extra_cap"], after["constraints"]["extra_cap"])
        for rid, rule in self.spec["reactions"].items():
            for field, value in reaction_fields(m.reactions.get_by_id(rid)).items():
                self.assertTrue(_same(field, value, rule["after"][field]), (rid, field))
        self.assertEqual(m.reactions.biomass_C.notes["canonical_trna_biomass_mode"], "split_v1")
        self.assertNotIn("metadata_previous_canonical_trna_biomass_mode", m.reactions.biomass_C.notes)
        self.assertNotIn("curated_gpr_correction", m.reactions.R612.notes)
        self.assertNotIn("gap_fill_direction_status", m.reactions.R_NTP1.notes)
        # SBML cannot carry arbitrary solver constraints; verify that constraint above,
        # then round-trip the representable model without claiming otherwise.
        m.remove_cons_vars(m.constraints.extra_cap)
        state = semantics(m)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "selected.xml"
            write_sbml_model(m, path)
            loaded = read_sbml_model(path)
        self.assertEqual(state, semantics(loaded))
        notes = annotations(loaded)
        self.assertTrue(all(r["status"] == "already_correct" for r in apply_metadata_reaction_selection(loaded)["items"]))
        self.assertEqual(notes, annotations(loaded))

    def test_second_application_does_not_accumulate_or_restore_old_notes(self):
        apply_metadata_reaction_selection(self.model)
        state = (semantics(self.model), annotations(self.model))
        report = apply_metadata_reaction_selection(self.model)
        self.assertTrue(all(r["status"] == "already_correct" for r in report["items"]))
        self.assertEqual(state, (semantics(self.model), annotations(self.model)))

    def test_metadata_selection_must_preserve_experimental_trna_biomass(self):
        m = self.model
        r = m.reactions.biomass_C
        before = (reaction_fields(r), dict(r.notes))
        self.assertEqual(self.spec["reactions"][r.id]["fields"], [])
        report = apply_metadata_reaction_selection(m)
        row = next(row for row in report["items"] if row["item"] == r.id)
        self.assertEqual(row["status"], "already_correct")
        self.assertEqual(before, (reaction_fields(r), r.notes))
        residues = {met: c for met, c in r.metabolites.items() if met.id.startswith("trna_biomass_residue_")}
        self.assertEqual(len(residues), 20)
        self.assertTrue(all(c < 0 for c in residues.values()))
        # The old mistake must fail the build status, even if the mode label survives.
        r.add_metabolites({met: 0 for met in residues}, combine=False)
        broken = reaction_fields(r)
        report = apply_metadata_reaction_selection(m)
        self.assertFalse(report["complete"])
        self.assertEqual(next(row["status"] for row in report["items"] if row["item"] == r.id), "conflict")
        self.assertEqual(broken, reaction_fields(r))

    def test_conflicting_equation_gpr_and_species_are_preserved_locally(self):
        m = self.model
        m.reactions.R159.gene_reaction_rule = "unknown_version_gene"
        m.reactions.R_NTP1.lower_bound = -20
        met = next(iter(m.reactions.biomass_C.metabolites))
        m.reactions.biomass_C.add_metabolites({met: 0.1})
        before = {rid: reaction_fields(m.reactions.get_by_id(rid)) for rid in ("R159", "R_NTP1", "biomass_C")}
        report = apply_metadata_reaction_selection(m)
        self.assertFalse(report["complete"])
        rows = {r["item"]: r for r in report["items"]}
        for rid, fields in before.items():
            self.assertEqual(rows[rid]["status"], "conflict")
            self.assertEqual(fields, reaction_fields(m.reactions.get_by_id(rid)))
        self.assertEqual(rows["R1302"]["status"], "applied")
        reaction = m.reactions.R2041
        fields = reaction_fields(reaction)
        next(iter(reaction.metabolites)).charge = 999
        report = apply_metadata_reaction_selection(m)
        row = next(row for row in report["items"] if row["item"] == reaction.id)
        self.assertEqual(row["status"], "conflict")
        self.assertEqual(fields, reaction_fields(reaction))

    def test_boolean_reordering_is_already_correct(self):
        r = self.model.reactions.R159
        target = self.spec["reactions"][r.id]["after"]["gpr"]
        r.gene_reaction_rule = " or ".join(reversed(target.split(" or ")))
        rule = r.gene_reaction_rule
        report = apply_metadata_reaction_selection(self.model)
        self.assertTrue(report["complete"])
        row = next(row for row in report["items"] if row["item"] == r.id)
        self.assertEqual(next(f["status"] for f in row["fields"] if f["field"] == "gpr"), "already_correct")
        self.assertEqual(rule, r.gene_reaction_rule)

    def test_exact_pre_export_states_without_general_numeric_tolerance(self):
        m = self.model
        rule = self.spec["reactions"]["R1372"]
        r = m.reactions.R1372
        r.add_metabolites({m.metabolites.get_by_id(mid): c for mid, c in rule["before_export"]["stoichiometry"].items()}, combine=False)
        for mid in self.spec["reactions"]["biomass_C"]["species_before_export"]:
            m.metabolites.get_by_id(mid).charge = None
        self.assertTrue(apply_metadata_reaction_selection(m)["complete"])
        met = next(iter(r.metabolites))
        r.add_metabolites({met: 1e-10})
        fields = reaction_fields(r)
        row = next(row for row in apply_metadata_reaction_selection(m)["items"] if row["item"] == r.id)
        self.assertEqual(row["status"], "conflict")
        self.assertEqual(fields, reaction_fields(r))


if __name__ == "__main__":
    unittest.main()
