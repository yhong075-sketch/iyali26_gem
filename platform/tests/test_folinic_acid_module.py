"""Folinic acid module: rejection, exact edits, rerun and SBML roundtrip, knockout logic, context, flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.folinic_acid_module import CURATION_PATH, NEW_ID, apply_folinic_acid_module
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_"
            "gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_"
            "20261103.xml")
SHMT_MI, SHMT_CY, GCVT, FAU1 = "YALI1D28569g", "YALI1E19557g", "YALI1F04082g", "YALI1D01489g"


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def groups(model):
    return {g.name: sorted(r.id for r in g.members) for g in model.groups}


def state(model):
    return semantics(model), annotations(model), gene_table(model), groups(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class FolinicAcidModuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_baseline()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        precondition = "precondition differs"
        cases = {
            "r748_gpr": ("R748", precondition), "r748_name": ("R748", precondition),
            "r748_annotation": ("R748", precondition), "r748_note": ("R748", precondition),
            "r990_bounds": ("R990", precondition), "r990_stoichiometry": ("R990", precondition),
            "species_charge": (None, "identity differs"), "m159_new_consumer": (None, "reaction context differs"),
            "gene_annotation": (None, "identity or reaction context"), "missing_group": ("R748", precondition),  # R748's pinned group membership fails first
            "thfat_exists_other": (None, "another definition"), "missing_reaction": ("R990", "R990 missing"),
        }
        baseline = state(self.model)
        for case, (_rid, message) in cases.items():
            bad = self.model.copy()
            if case == "r748_gpr":
                bad.reactions.R748.gene_reaction_rule = SHMT_MI
            elif case == "r748_name":
                bad.reactions.R748.name = "methenyl-THF hydrolase"
            elif case == "r748_annotation":
                bad.reactions.R748.annotation = {**bad.reactions.R748.annotation, "ec-code": "2.1.2.1"}
            elif case == "r748_note":
                bad.reactions.R748.notes = {**bad.reactions.R748.notes, "unreviewed_note": "added elsewhere"}
            elif case == "r990_bounds":
                bad.reactions.R990.bounds = (-1000.0, 1000.0)
            elif case == "r990_stoichiometry":
                bad.reactions.R990.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif case == "species_charge":
                bad.metabolites.get_by_id("m160[C_cy]").charge = -1
            elif case == "m159_new_consumer":
                bad.reactions.R1907.add_metabolites({bad.metabolites.get_by_id("m159[C_cy]"): 1})
            elif case == "gene_annotation":
                # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
                gene = bad.genes.get_by_id(SHMT_CY)
                gene.annotation = {**gene.annotation, "ncbigene": "2900000"}
            elif case == "missing_group":
                bad.remove_groups([g for g in bad.groups if g.name == "One carbon pool by folate"])
            elif case == "thfat_exists_other":
                other = bad.reactions.R748.copy()
                other.id = NEW_ID
                bad.add_reactions([other])
            else:
                bad.remove_reactions([bad.reactions.R990])
            before = state(bad)
            with self.subTest(case=case), self.assertRaisesRegex(ValueError, message):
                apply_folinic_acid_module(bad)
            self.assertEqual(state(bad), before, case)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_partial_state_is_rejected(self):
        bad = self.model.copy()
        bad.reactions.R990.bounds = (0.0, 0.0)
        bad.reactions.R990.notes = dict(self.spec["reactions"]["R990"]["after"]["notes"])
        before = state(bad)
        with self.assertRaisesRegex(ValueError, "Only part"):
            apply_folinic_acid_module(bad)
        self.assertEqual(state(bad), before)

    def test_unexpected_scope_is_rejected(self):
        edits = {
            "status": lambda s: s.update(status="draft"),
            "curation_id": lambda s: s.update(curation_id="FOLINIC-ACID-MODULE"),
            "r748_stoichiometry": lambda s: s["reactions"]["R748"]["after"]["stoichiometry"].update({"m28[C_mi]": 2.0}),
            "r990_open": lambda s: s["reactions"]["R990"]["after"].update(bounds=[0.0, 1000.0]),
            "r990_gpr": lambda s: s["reactions"]["R990"]["after"].update(gpr=FAU1),
            "new_rule": lambda s: s["new_reaction"].update(gpr=SHMT_CY),
            "new_reversible": lambda s: s["new_reaction"].update(bounds=[-1000.0, 1000.0]),
            "extra_gene": lambda s: s["genes"].update({FAU1: copy.deepcopy(s["genes"][SHMT_CY])}),
        }
        for name, edit in edits.items():
            spec = copy.deepcopy(self.spec)
            edit(spec)
            with self.subTest(edit=name), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_folinic_acid_module(self.model.copy(), spec)

    def test_changes_only_the_targets_and_survives_roundtrip(self):
        model = read_baseline()
        before_ids = [r.id for r in model.reactions]
        before_state = state(model)
        result = apply_folinic_acid_module(model)
        self.assertEqual(result["status"], "applied")
        after_state = state(model)
        # Only R748, R990 and the new R_THFAT differ in reactions, annotations and notes.
        sem_b, ann_b, genes_b, groups_b = before_state
        sem_a, ann_a, genes_a, groups_a = after_state
        self.assertEqual(sem_a["species"], sem_b["species"])
        self.assertEqual(sem_a["genes"], sem_b["genes"], "No gene may be added or removed")
        changed = {r for r in set(sem_a["reactions"]) | set(sem_b["reactions"])
                   if sem_a["reactions"].get(r) != sem_b["reactions"].get(r) or ann_a.get(r) != ann_b.get(r)}
        self.assertEqual(changed, {"R748", "R990", NEW_ID})
        self.assertEqual([r.id for r in model.reactions][:-1], before_ids)
        self.assertEqual({g for g in genes_a if genes_a[g] != genes_b[g]}, {SHMT_MI, SHMT_CY, GCVT})
        self.assertEqual({g for g in groups_a if groups_a[g] != groups_b[g]}, {"One carbon pool by folate"})
        self.assertEqual(sorted(set(groups_a["One carbon pool by folate"]) - set(groups_b["One carbon pool by folate"])),
                         [NEW_ID])
        new = model.reactions.get_by_id(NEW_ID)
        self.assertEqual(new.check_mass_balance(), model.reactions.R748.check_mass_balance())
        self.assertEqual(new.check_mass_balance(), {"charge": -1.0})
        self.assertEqual(apply_folinic_acid_module(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "folinic.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertEqual(state(loaded)[1:], after_state[1:])
        self.assertEqual(apply_folinic_acid_module(loaded)["status"], "already_correct")

    def test_knockout_logic(self):
        model = self.model.copy()
        self.assertEqual(model.reactions.R748.gene_reaction_rule, GCVT)
        apply_folinic_acid_module(model)
        cases = {(GCVT,): False, (SHMT_MI,): False, (SHMT_CY,): False, (SHMT_MI, SHMT_CY): True}
        for genes, closes in cases.items():
            with model, self.subTest(genes=genes):
                for gene_id in genes:
                    model.genes.get_by_id(gene_id).knock_out()
                for rid in ("R748", NEW_ID, "R358", "R359"):
                    self.assertEqual(model.reactions.get_by_id(rid).bounds == (0, 0), closes, rid)
        # R990 is closed and gene-less; R64 stays on the FAU1-type gene alone.
        self.assertEqual(model.reactions.R990.bounds, (0.0, 0.0))
        self.assertEqual(model.reactions.R64.gene_reaction_rule, FAU1)

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_folinic_acid_module(model)["status"], "applied")
            self.assertIn(NEW_ID, model.reactions)
        self.assertNotIn(NEW_ID, model.reactions)
        self.assertEqual(state(model), before, "Context exit left curation edits behind")
        self.assertEqual(apply_folinic_acid_module(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).folinic_acid_module)
        self.assertTrue(parse_args(["--folinic-acid-module"]).folinic_acid_module)


if __name__ == "__main__":
    unittest.main()
