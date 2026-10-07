"""The opt-in hypothesis requires three genes without editing reaction chemistry."""
import copy
import itertools
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Model, Metabolite, Reaction
from cobra.core.gene import GPR
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.patches import VATPASE_HYPOTHESIS_PATH, apply_vatpase_gpr_hypothesis
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics


class VatpaseHypothesisTests(unittest.TestCase):
    def setUp(self):
        self.spec = json.loads(VATPASE_HYPOTHESIS_PATH.read_text())
        self.model = Model('vatpase_hypothesis')
        for rid, row in self.spec['reactions'].items():
            for mid, attributes in row['species'].items():
                if mid not in self.model.metabolites:
                    self.model.add_metabolites([Metabolite(mid, name=mid, **attributes)])
            reaction = Reaction(rid)
            reaction.add_metabolites({self.model.metabolites.get_by_id(mid): c for mid, c in row['stoichiometry'].items()})
            reaction.bounds = row['bounds']
            reaction.gene_reaction_rule = row['before_gpr']
            reaction.notes = {'existing_note': 'preserve'}
            self.model.add_reactions([reaction])
        for gene in self.model.genes:
            gene.name = gene.id
        self.model.objective = self.model.reactions.R794

    def test_exact_boolean_requirement_and_roundtrip(self):
        before = semantics(self.model)
        self.assertEqual([r['status'] for r in apply_vatpase_gpr_hypothesis(self.model)], ['applied'] * 2)
        required = set(self.spec['required_genes'])
        for rid, row in self.spec['reactions'].items():
            reaction = self.model.reactions.get_by_id(rid)
            old = GPR.from_string(row['before_gpr'])
            genes = sorted(old.genes)
            for bits in itertools.product((False, True), repeat=len(genes)):
                knocked = {g for g, absent in zip(genes, bits) if absent}
                self.assertEqual(reaction.gpr.eval(knocked), not (required & knocked) and old.eval(knocked))
            for gid in required:
                self.assertFalse(reaction.gpr.eval({gid}))
            self.assertEqual(reaction.notes['existing_note'], 'preserve')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'candidate.xml'
            write_deterministic_sbml_model(self.model, path)
            loaded = read_sbml_model(str(path))
        self.assertEqual([r['status'] for r in apply_vatpase_gpr_hypothesis(loaded)], ['already_correct'] * 2)
        for rid, row in self.spec['reactions'].items():
            self.assertEqual(loaded.reactions.get_by_id(rid).notes, self.model.reactions.get_by_id(rid).notes)
            loaded.reactions.get_by_id(rid).gene_reaction_rule = row['before_gpr']
        self.assertEqual(semantics(loaded), before)

    def test_second_reaction_conflict_does_not_partially_edit_first(self):
        self.model.reactions.R795.upper_bound = 1
        before = semantics(self.model)
        notes = copy.deepcopy(self.model.reactions.R794.notes)
        with self.assertRaisesRegex(ValueError, 'precondition differs for R795'):
            apply_vatpase_gpr_hypothesis(self.model)
        self.assertEqual(semantics(self.model), before)
        self.assertEqual(self.model.reactions.R794.notes, notes)
        self.model.reactions.R795.bounds = self.spec['reactions']['R795']['bounds']
        self.model.reactions.R794.gene_reaction_rule = self.spec['reactions']['R794']['after_gpr']
        with self.assertRaisesRegex(ValueError, 'Mixed V-ATPase hypothesis state'):
            apply_vatpase_gpr_hypothesis(self.model)
        self.assertEqual(self.model.reactions.R795.gene_reaction_rule, self.spec['reactions']['R795']['before_gpr'])

    def test_opt_in_and_separate_output_guards(self):
        self.assertFalse(parse_args([]).vatpase_gpr_hypothesis)
        self.assertTrue(parse_args(['--vatpase-gpr-hypothesis']).vatpase_gpr_hypothesis)
        with self.assertRaisesRegex(ValueError, 'requires an offline/no-solve metadata build'):
            build_reference_chain(vatpase_gpr_hypothesis=True)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'existing.xml'
            path.write_text('preserved')
            with self.assertRaises(FileExistsError):
                build_reference_chain(vatpase_gpr_hypothesis=True, no_solve=True, allow_network=False, output_model_path=path)
            self.assertEqual(path.read_text(), 'preserved')


if __name__ == '__main__':
    unittest.main()
