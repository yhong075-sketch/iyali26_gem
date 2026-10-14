"""In-memory checks for the 2026-10-13 R132 removal; no optimization."""

import copy
import json
import unittest

from cobra import Gene, Metabolite, Model, Reaction

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.r132_mito_adh import SPEC_PATH, apply_r132_mito_adh

SPEC = json.loads(SPEC_PATH.read_text())
# Other reactions that use R132 metabolites or genes in the baseline, reduced to what the guard reads.
OTHER_USERS = {"R135": ["m278[C_mi]"], "R2070": ["m278[C_mi]"], "R130": [], "R257": []}
OTHER_GENES = {"YALI1A15227g": "R130", "YALI1E15452g": "R257"}


def _reaction(model, rid, stoichiometry, bounds, gpr):
    reaction = Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([reaction])
    for mid in stoichiometry:
        if mid not in model.metabolites:
            model.add_metabolites([Metabolite(mid, compartment=mid.split("[")[1][:-1])])
    reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in stoichiometry.items()})
    reaction.gene_reaction_rule = gpr
    return reaction


def fixture():
    model = Model("r132_fixture")
    _reaction(model, "R132", SPEC["before"]["stoichiometry"], SPEC["before"]["bounds"], SPEC["before"]["gpr"])
    for rid, state in SPEC["kept_reactions"].items():
        _reaction(model, rid, state["stoichiometry"], state["bounds"], state["gpr"])
    for rid, mids in OTHER_USERS.items():
        gene = next((g for g, r in OTHER_GENES.items() if r == rid), "")
        _reaction(model, rid, {mid: -1.0 for mid in mids}, (0.0, 1000.0), gene)
    for gid, row in SPEC["genes"].items():
        if row["model_refseq"]:
            model.genes.get_by_id(gid).annotation["refseq"] = row["model_refseq"]
    return model


class R132MitoAdhTests(unittest.TestCase):
    def test_removes_r132_only_and_is_idempotent(self):
        model = fixture()
        genes, metabolites = {g.id for g in model.genes}, {m.id for m in model.metabolites}
        kept = {rid: model_definition(model)["reactions"][rid] for rid in SPEC["kept_reactions"]}
        result = apply_r132_mito_adh(model)
        self.assertEqual(result["status"], "applied")
        self.assertNotIn("R132", model.reactions)
        self.assertEqual({g.id for g in model.genes}, genes)  # Genes are kept.
        self.assertEqual({m.id for m in model.metabolites}, metabolites)  # Metabolites are kept.
        self.assertEqual(result["orphaned_genes"],
                         sorted(g for g, row in SPEC["genes"].items() if row["only_in_r132"]))
        self.assertEqual(result["metabolite_usage_after"], {"m279[C_mi]": ["R1114"]})
        self.assertEqual({rid: model_definition(model)["reactions"][rid] for rid in kept}, kept)
        definition = model_definition(model)
        self.assertEqual(apply_r132_mito_adh(model)["status"], "already_correct")
        self.assertEqual(model_definition(model), definition)

    def test_stale_inputs_fail_closed_without_edits(self):
        def bounds(m): m.reactions.R132.bounds = (0.0, 1000.0)
        def gpr(m): m.reactions.R132.gene_reaction_rule = "YALI1A15227g"
        def kept(m): m.reactions.R1114.bounds = (0.0, 1000.0)
        def refseq(m): m.genes.YALI1F11481g.annotation["refseq"] = "XP_505149.1"
        def ethanol_user(m): _reaction(m, "extra", {"m279[C_mi]": -1.0}, (0.0, 1000.0), "")
        def membership(m): _reaction(m, "extra", {}, (0.0, 1000.0), "YALI1F11481g")
        for change in (bounds, gpr, kept, refseq, ethanol_user, membership):
            with self.subTest(change.__name__):
                model = fixture()
                change(model)
                before = model_definition(model)
                with self.assertRaises(ValueError):
                    apply_r132_mito_adh(model)
                self.assertEqual(model_definition(model), before)

    def test_unexpected_scope_is_rejected(self):
        spec = copy.deepcopy(SPEC)
        spec["reaction"] = "R130"
        with self.assertRaises(ValueError):
            apply_r132_mito_adh(fixture(), spec)

    def test_flag_requires_the_respiratory_package(self):
        self.assertFalse(parse_args([]).r132_mito_adh)
        self.assertTrue(parse_args(["--r132-mito-adh"]).r132_mito_adh)
        with self.assertRaises(ValueError):
            build_reference_chain(r132_mito_adh=True)


if __name__ == "__main__":
    unittest.main()
