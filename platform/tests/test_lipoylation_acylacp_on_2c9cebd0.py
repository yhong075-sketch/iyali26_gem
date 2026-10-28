"""In-memory checks of the lipoylation and acyl-ACP steps on the main candidate 2c9cebd0 (2026-10-23); no optimization."""

import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.acyl_acp_chemistry import SPEC_PATH as ACYL_SPEC_PATH, apply_acyl_acp_chemistry
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.heme_iron_biomass import apply_heme_iron_biomass
from scripts.gem_annotate.lipoylation_module import SPEC_PATH as LIPO_SPEC_PATH, apply_lipoylation_module
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model

LIPO = json.loads(LIPO_SPEC_PATH.read_text())
ACYL = json.loads(ACYL_SPEC_PATH.read_text())
BASELINE = "E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml"


def snapshot(model):
    return (model_definition(model), sorted(g.id for g in model.genes),
            {r.id: dict(r.notes) for r in model.reactions})


class LipoylationAcylAcpOnMainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.baseline = read_sbml_model(str(MODEL.candidate_file(BASELINE)))

    def test_both_steps_change_only_their_targets(self):
        model = self.baseline.copy()
        before = model_definition(model)
        self.assertEqual(apply_lipoylation_module(model)["status"], "applied")
        self.assertEqual(apply_acyl_acp_chemistry(model)["status"], "applied")
        after = model_definition(model)
        # Union of keys, so a removed reaction or species also counts as a change.
        changed_reactions = sorted(r for r in set(before["reactions"]) | set(after["reactions"])
                                   if before["reactions"].get(r) != after["reactions"].get(r))
        self.assertEqual(changed_reactions, sorted({"R362", "R364", "biomass_C", *LIPO["new_reactions"],
                                                    *ACYL["changed_reactions"]}))
        changed_species = sorted(m for m in set(before["metabolites"]) | set(after["metabolites"])
                                 if before["metabolites"].get(m) != after["metabolites"].get(m))
        self.assertEqual(changed_species, sorted({*LIPO["formula_fills"], *LIPO["new_metabolites"],
                                                  *ACYL["formula_fills"]}))
        before_genes, after_genes = {g.id for g in self.baseline.genes}, {g.id for g in model.genes}
        self.assertEqual(before_genes - after_genes, set())
        self.assertEqual(sorted(after_genes - before_genes), sorted(LIPO["new_genes"]))
        self.assertEqual(after["objective"], before["objective"])

    def test_reruns_are_no_ops_after_sbml_round_trip(self):
        model = self.baseline.copy()
        apply_lipoylation_module(model)
        apply_acyl_acp_chemistry(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "combined.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        state = snapshot(loaded)
        self.assertEqual(apply_lipoylation_module(loaded)["status"], "already_correct")
        self.assertEqual(apply_acyl_acp_chemistry(loaded)["status"], "already_correct")
        self.assertEqual(snapshot(loaded), state)

    def test_heme_guard_rejects_lipoylated_biomass(self):
        # Why the builder runs the lipoylation step after the heme/iron step: the heme guard pins biomass_C.
        model = self.baseline.copy()
        self.assertEqual(apply_heme_iron_biomass(model)["status"], "already_correct")
        apply_lipoylation_module(model)
        with self.assertRaisesRegex(ValueError, "biomass_C composition differs from the reviewed baseline"):
            apply_heme_iron_biomass(model)


if __name__ == "__main__":
    unittest.main()
