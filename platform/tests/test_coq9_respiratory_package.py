"""In-memory checks for the 2026-10-08 CoQ9/respiratory package; no optimization."""

import copy
import json
import unittest

from cobra import Gene, Metabolite, Model, Reaction

from scripts.gem_annotate.coq9 import FUNCTIONAL_GPR_PATH, exact_residual
from scripts.gem_annotate.coq9_respiratory_package import PACKAGE_PATH, apply_coq9_respiratory_package
from scripts.gem_annotate.energy_candidates import model_definition

SPEC = json.loads(PACKAGE_PATH.read_text())
R695 = json.loads(FUNCTIONAL_GPR_PATH.read_text())
COUNTERPARTS = {
    "R247": ({"m10[C_cy]": -1.0, "m445[C_cy]": -1.0, "m446[C_cy]": -1.0, "m202[C_cy]": 1.0, "m443[C_cy]": 1.0},
             (0.0, 1000.0), "YALI1E10364g or YALI1F18920g or YALI1D33380g"),
    "R348": ({"m10[C_cy]": -2.0, "m123[C_cy]": -1.0, "m456[C_cy]": -1.0, "m122[C_cy]": 1.0, "m562[C_cy]": 1.0},
             (0.0, 1000.0), "YALI1B04433g"),
}


def _metabolite(model, mid, identity=None):
    if mid in model.metabolites:
        return model.metabolites.get_by_id(mid)
    identity = identity or {"formula": None, "charge": None, "compartment": mid.split("[")[1][:-1]}
    met = Metabolite(mid, **identity)
    model.add_metabolites([met])
    return met


def _reaction(model, rid, stoichiometry, bounds, gpr, species=None):
    reaction = Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([reaction])
    reaction.add_metabolites({_metabolite(model, mid, (species or {}).get(mid)): c
                              for mid, c in stoichiometry.items()})
    reaction.gene_reaction_rule = gpr
    return reaction


def fixture():
    model = Model("package_fixture")
    species = {}
    for item in SPEC["items"]:
        for row in item.get("set", {}).values():
            species.update(row["after_species"])
    for mid, identity in species.items():
        _metabolite(model, mid, identity)
    for item in SPEC["items"]:
        for rid, row in item.get("set", {}).items():
            _reaction(model, rid, row["before"]["stoichiometry"], row["before"]["bounds"], row["before"]["gpr"])
        for rid, state in item.get("remove", {}).items():
            _reaction(model, rid, state["stoichiometry"], state["bounds"], state["gpr"])
        for gid, identity in item.get("genes", {}).items():
            if gid not in model.genes:  # Orphan genes in the real model (e.g. COQ6/COQ8 candidates).
                model.genes.append(Gene(gid))
            model.genes.get_by_id(gid).annotation.update(identity)
    guard = R695["guard"]
    _reaction(model, "R695", {mid: s["coefficient"] for mid, s in guard["species"].items()},
              guard["bounds"], guard["gpr"],
              {mid: {k: s[k] for k in ("formula", "charge", "compartment")} for mid, s in guard["species"].items()})
    model.genes.append(Gene("YALI1F34675g"))
    for gid, identity in R695["genes"].items():
        model.genes.get_by_id(gid).annotation["refseq"] = identity["refseq"]
    for rid, (stoichiometry, bounds, gpr) in COUNTERPARTS.items():
        _reaction(model, rid, stoichiometry, bounds, gpr)
    for item in SPEC["items"]:
        for rid, row in item.get("set_notes", {}).items():
            if rid not in model.reactions:
                _reaction(model, rid, {}, (0.0, 1000.0), "")
            model.reactions.get_by_id(rid).notes[row["key"]] = row["before"]
    model.objective = "R695"
    return model


class CoQ9RespiratoryPackageTests(unittest.TestCase):
    def test_applies_curated_targets_and_is_idempotent(self):
        model = fixture()
        genes_before = {g.id for g in model.genes}
        result = apply_coq9_respiratory_package(model)
        self.assertEqual([r["status"] for r in result["items"]], ["applied"] * 7)
        for item in SPEC["items"]:
            for rid, row in item.get("set", {}).items():
                reaction = model.reactions.get_by_id(rid)
                self.assertEqual({m.id: c for m, c in reaction.metabolites.items()}, row["after"]["stoichiometry"])
                self.assertEqual(list(reaction.bounds), row["after"]["bounds"])
                self.assertEqual(exact_residual(reaction), {})
            for rid in item.get("remove", {}):
                self.assertNotIn(rid, model.reactions)
        self.assertNotIn("m108[C_cy]", model.metabolites)
        self.assertNotIn("m110[C_cy]", model.metabolites)
        self.assertEqual({m.compartment for m in model.reactions.R39.metabolites}, {"C_mi"})
        self.assertEqual(model.reactions.R695.gene_reaction_rule, "YALI1E18269g and YALI1F34675g")
        self.assertNotIn("YALI1F32476g", {g.id for g in model.reactions.R2062.genes})
        self.assertNotIn("YALI1E06573g", {g.id for g in model.reactions.R2062.genes})
        b1 = next(item for item in SPEC["items"] if item["id"] == "B1")
        self.assertEqual(model.reactions.R570.notes["complex_i_scope"], b1["set_notes"]["R570"]["after"])
        self.assertEqual({g.id for g in model.genes}, genes_before)  # Genes are never removed.
        definition = model_definition(model)
        again = apply_coq9_respiratory_package(model)
        self.assertEqual({r["status"] for r in again["items"]}, {"already_correct"})
        self.assertEqual(model_definition(model), definition)

    def test_stale_target_fails_closed_without_edits(self):
        model = fixture()
        model.reactions.R304.upper_bound = 999.0
        before = model_definition(model)
        with self.assertRaises(ValueError):
            apply_coq9_respiratory_package(model)
        self.assertEqual(model_definition(model), before)

    def test_stale_r695_fails_before_any_edit(self):
        model = fixture()
        model.reactions.R695.gene_reaction_rule = "YALI1F34675g"
        before = model_definition(model)
        with self.assertRaises(ValueError):
            apply_coq9_respiratory_package(model)
        self.assertEqual(model_definition(model), before)

    def test_unexpected_user_of_orphan_metabolite_fails_before_any_edit(self):
        model = fixture()
        _reaction(model, "extra_m108_user", {"m108[C_cy]": -1.0}, (0.0, 1000.0), "")
        before = model_definition(model)
        with self.assertRaises(ValueError):
            apply_coq9_respiratory_package(model)
        self.assertEqual(model_definition(model), before)

    def test_duplicate_removal_requires_exact_counterpart(self):
        model = fixture()
        model.reactions.R247.add_metabolites({model.metabolites.get_by_id("m10[C_cy]"): 1.0})
        before = model_definition(model)
        with self.assertRaises(ValueError):
            apply_coq9_respiratory_package(model)
        self.assertEqual(model_definition(model), before)

    def test_unexpected_scope_is_rejected(self):
        spec = copy.deepcopy(SPEC)
        spec["items"] = spec["items"][:-1]
        with self.assertRaises(ValueError):
            apply_coq9_respiratory_package(fixture(), spec)


if __name__ == "__main__":
    unittest.main()
