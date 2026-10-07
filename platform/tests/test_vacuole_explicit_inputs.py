"""Explicit input failures and archived-run audits; zero optimization calls."""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from cobra import Model, Reaction, Metabolite
from scripts.gem_annotate.execution import execution_limits
from tools.validate_vacuole_supply import (
    ROOT, audit_saved_run, bind_inputs, input_identity, main, math_signature, sha,
)

from scripts.gem_annotate.config import load_project_paths
from scripts.gem_annotate.model_layout import SCRATCH_DIR
# Both saved task folders live in the research workspace, not in git.
TASK = load_project_paths().task_outputs / 'dipeptide_chemistry_routes_20260924'
PREVIOUS = load_project_paths().task_outputs / 'vacuole_open_supply_20260924'


class ExplicitInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=SCRATCH_DIR)
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.base, self.candidate = self.folder/'base.xml', self.folder/'candidate.xml'
        self.base.write_text('baseline')
        self.candidate.write_text('candidate')
        self.guard = execution_limits(no_solve=True, allow_network=False)
        self.attempts = self.guard.__enter__()
        self.addCleanup(self.guard.__exit__, None, None, None)

    def tearDown(self):
        self.assertEqual(self.attempts, {'optimization': 0, 'network': 0})

    def test_explicit_paths_sha_and_manifest_binding(self):
        good = bind_inputs(self.base, sha(self.base), self.candidate, sha(self.candidate))
        self.assertEqual(good['candidate']['path'], str(self.candidate))
        for path, digest, kind in [(self.folder/'missing.xml', sha(self.candidate), FileNotFoundError),
                                   (self.candidate, '0'*64, ValueError),
                                   (self.candidate, None, ValueError)]:
            with self.assertRaises(kind):
                input_identity(path, digest)
        manifest = self.folder/'candidate.build.json'
        record = {'source': str(self.base), 'source_sha256': sha(self.base),
                  'output': str(self.candidate), 'output_sha256': sha(self.candidate)}
        manifest.write_text(json.dumps(record))
        self.assertTrue(bind_inputs(self.base, sha(self.base), self.candidate,
            build_manifest=manifest)['build_manifest']['source_and_output_bind_requested_paths'])
        # Identical bytes at another path are still not the built file requested by this manifest.
        other = self.folder/'other.xml'
        other.write_bytes(self.candidate.read_bytes())
        with self.assertRaisesRegex(ValueError, 'different candidate path'):
            bind_inputs(self.base, sha(self.base), other, build_manifest=manifest)
        record['source_sha256'] = '0'*64
        manifest.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, 'different baseline'):
            bind_inputs(self.base, sha(self.base), self.candidate, build_manifest=manifest)

    def test_cli_never_falls_back_to_config_candidate(self):
        args = ['--mode', 'fresh_validation', '--config', str(PREVIOUS/'config.json'),
                '--baseline-model', str(self.base), '--baseline-sha256', sha(self.base),
                '--output', str(self.folder/'output'), '--budget', str(self.folder/'budget.json'),
                '--expected-diff', 'four_connections']
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(args)
        with self.assertRaises(FileNotFoundError):
            main(args + ['--candidate-model', str(self.folder/'missing.xml'),
                         '--candidate-sha256', '0'*64])
        self.assertFalse((self.folder/'output').exists())
        self.assertFalse((self.folder/'budget.json').exists())

    def test_mathematical_signature_ignores_chemistry_but_detects_constraints(self):
        model = Model('test')
        metabolite = Metabolite('m', name='missing formula', compartment='c')
        reaction = Reaction('r', lower_bound=0, upper_bound=1)
        reaction.add_metabolites({metabolite: 1})
        model.add_reactions([reaction])
        before = math_signature(model)
        metabolite.name, metabolite.formula, metabolite.charge = 'glycine', 'C2H5NO2', 0
        self.assertEqual(math_signature(model), before)
        model.add_cons_vars(model.problem.Constraint(reaction.flux_expression, ub=0.5, name='extra'))
        self.assertNotEqual(math_signature(model), before)

    @unittest.skipUnless((PREVIOUS / 'run').is_dir() and (TASK / 'previous_source_snapshot').is_dir(),
                         'saved task folders not found; set IYALI26_RESEARCH_ROOT')
    def test_audit_saved_run_uses_archived_executed_source_without_optimization(self):
        config = json.loads((PREVIOUS/'config.json').read_text())
        identities = bind_inputs(config['model'], config['model_sha256'], config['candidate'],
                                  build_manifest=PREVIOUS/'E5_vacuole_open.build.json')
        with patch('cobra.Model.optimize', side_effect=AssertionError('audit must not optimize')):
            result = audit_saved_run(PREVIOUS/'run', TASK/'previous_source_snapshot',
                self.folder/'audit', PREVIOUS/'budget.json', PREVIOUS/'config.json', identities)
        self.assertEqual(result['optimization_calls_this_audit'], 0)
        self.assertEqual(result['saved_optimization_calls'], 114)
        source = next(r for r in result['source_identity'] if r['path'] == 'scripts/validate_vacuole_supply.py')
        self.assertFalse(source['current_matches_executed'])
        self.assertEqual(source['executed_sha256'], source['archived_sha256'])


if __name__ == '__main__':
    unittest.main()
