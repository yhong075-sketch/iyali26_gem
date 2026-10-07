"""The opt-in pipeline keeps the reviewed chemistry and fails on real drift."""
import copy
import json
import tempfile
import unittest
from html import unescape
from pathlib import Path

from cobra import Gene, Metabolite, Model, Reaction

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.coq9 import LITERATURE_PATH, apply_coq_literature_revision
from scripts.gem_annotate.execution import execution_limits
from scripts.gem_annotate.main import build_reference_chain
from tests.test_coq9_curation import annotations, semantics


class CoQLiteraturePipelineTest(unittest.TestCase):
    def test_opt_in_and_execution_guards(self):
        self.assertFalse(parse_args([]).coq_literature_revision)
        self.assertTrue(parse_args(["--coq-literature-revision"]).coq_literature_revision)
        with tempfile.TemporaryDirectory() as directory:
            options = dict(coq_literature_revision=True, no_solve=True,
                           allow_network=False, output_model_path=Path(directory) / "new.xml")
            for conflict in ({"no_solve": False}, {"allow_network": True},
                             {"coq9_mode": "qcycle"}, {"coq9_mode": "off"},
                             {"energy_candidate": "E5"}, {"vatpase_gpr_hypothesis": True},
                             {"provisional_capacity_path": Path("unused")},
                             {"trna_biomass_mode": "split"}, {"r608_curation_path": Path("unused")}):
                with self.subTest(conflict=conflict), self.assertRaises(ValueError):
                    build_reference_chain(**{**options, **conflict})
            options["output_model_path"].write_text("preserve")
            with self.assertRaises(FileExistsError):
                build_reference_chain(**options)
            self.assertEqual(options["output_model_path"].read_text(), "preserve")

    def test_serialization_equivalence_and_atomic_rejection(self):
        spec = json.loads(LITERATURE_PATH.read_text())
        with execution_limits(no_solve=True, allow_network=False) as attempts:
            model = Model("coq_context")
            rows = {**{rid: row["before"] for rid, row in spec["reactions"].items()},
                    **spec["context_guards"]}
            for rid, row in rows.items():
                for mid, props in row["species"].items():
                    if mid not in model.metabolites:
                        model.add_metabolites([Metabolite(mid, formula=props["formula"],
                            charge=props["charge"], compartment=props["compartment"])])
                reaction = Reaction(rid, name=row["name"])
                reaction.bounds = row["bounds"]
                reaction.add_metabolites({model.metabolites.get_by_id(mid): p["coefficient"]
                                         for mid, p in row["species"].items()})
                reaction.gene_reaction_rule = row["gpr"]
                reaction.annotation = copy.deepcopy(row["annotation"])
                reaction.notes = copy.deepcopy(row["notes"])
                model.add_reactions([reaction])
            # A species used only after the correction must already exist upstream.
            for row in spec["reactions"].values():
                for mid, props in row["after"]["species"].items():
                    if mid not in model.metabolites and mid not in spec["new_metabolites"]:
                        model.add_metabolites([Metabolite(mid, formula=props["formula"],
                            charge=props["charge"], compartment=props["compartment"])])
            for gid, identity in spec["genes"].items():
                if gid not in model.genes:
                    model.genes.append(Gene(gid))
                model.genes.get_by_id(gid).annotation["refseq"] = identity["refseq"]
            state = semantics(model), annotations(model), copy.deepcopy(model.notes)
            self.assertEqual(apply_coq_literature_revision(model)["status"], "disabled")
            self.assertEqual((semantics(model), annotations(model), model.notes), state)
            for kind in ("ec", "note", "chemistry"):
                bad = model.copy()
                if kind == "ec":
                    bad.reactions.R695.annotation["ec-code"] = "wrong"
                elif kind == "note":
                    bad.reactions.R385.notes["terminal_redox_convention"] = "unverified replacement"
                else:
                    bad.reactions.R18.upper_bound = 5
                before = semantics(bad), annotations(bad), copy.deepcopy(bad.notes)
                with self.subTest(kind=kind), self.assertRaises(ValueError):
                    apply_coq_literature_revision(bad, True)
                self.assertEqual((semantics(bad), annotations(bad), bad.notes), before)
            # In-memory pipeline values differ from SBML-read values only in encoding.
            for reaction in model.reactions:
                reaction.annotation = {k: [v] if isinstance(v, str) else v
                                       for k, v in reaction.annotation.items()}
                reaction.notes = {k: unescape(v) for k, v in reaction.notes.items()}
            self.assertEqual(apply_coq_literature_revision(model, True)["status"], "applied")
            self.assertEqual(apply_coq_literature_revision(model, True)["status"], "already_correct")
            self.assertEqual(attempts, {"optimization": 0, "network": 0})


if __name__ == "__main__":
    unittest.main()
