# E5 ATP energy repair candidate

The delivered model is [candidates/E5.xml](candidates/E5.xml). It remains an
explicit candidate; this delivery does not replace `model.xml` or activate E5
in the default build. The original input is preserved as
`model_metadata_trna_r1159_leak.xml` at the repository root.

Eight reaction definitions are repaired through
`data/energy_candidate_repairs.json`: R_PGAM1_PhosHydro, R_NTP3pp, R_NTP7, r0242,
R_NDP1, R72, R_CAT2p and R_OAADCm. Exact reaction/species checks reject conflicting
inputs before editing. Candidate definitions persist through metadata selection
and SBML export/reload. GPRs, medium and biomass requirements are unchanged by
this energy patch.

## Observed results

Under the recorded closed-system protocol, the maximum ATP, GTP, UTP and CTP
dissipation was zero, with `optimal` status for all four objectives and a
1e-7 tolerance. The exact zero-flux case was feasible. Under the original
PO1f/SD-Leu/plasmid-selection context, biomass flux was 1.429177892279185 h^-1
versus 1.8718823069403008 h^-1 for E0 (23.6502% lower). Normal ATP synthesis
remained active (R171 flux 53.485020907031746), and the maintenance lower bound
remained 7.8625.

These observations support removal of the tested energy-regeneration routes.
They do not establish universal thermodynamic consistency, native GPR validity,
compartment identity, or agreement with experimental growth. In particular,
R_CAT2p localization and other unresolved chemistry remain limitations.

The original repair run recorded 58 optimization calls and 23 passing tests.
Its [final E5 solve records](E5_validation/manifest.json),
[energy results](E5_validation/energy_iterations.tsv),
[growth comparison](growth_comparison.tsv),
[reaction repairs](accepted_candidate_fixes.tsv),
[chemical review](chemistry/REPORT.md) and [audit](AUDIT.md) are included.
Historical records retain their original local paths, code hashes and broader
working-tree context. This delivery includes the final E5 full solve records
and selected historical evidence, not every intermediate artifact referenced
by those records. The E3 XML is included as the known-anomaly regression control.

## Rebuild and check from a clean checkout

Install the dependencies in `pyproject.toml`; optimization requires a working
Gurobi license. Run from the repository root. Use new output paths on each run.

```bash
python -B scripts/build_energy_candidates.py \
  --variants E5 --output-dir artifacts/energy_rebuild

IYALI26_ENERGY_BUDGET=artifacts/atp_candidate_repair_20260924/new_test_budget.json \
IYALI26_ENERGY_TEST_OUTPUT=artifacts/atp_candidate_repair_20260924/new_behavior_tests \
python -B -m unittest tests.test_energy_candidate_build tests.test_energy_candidate_behavior

python -B -m scripts.validate_energy_candidates \
  --model E5=artifacts/atp_candidate_repair_20260924/candidates/E5.xml \
  --other-energy E5 \
  --budget artifacts/atp_candidate_repair_20260924/new_validation_budget.json \
  --output artifacts/atp_candidate_repair_20260924/new_validation
```

The standalone builder starts from the pinned completed input and performs no
optimization or network requests. It verifies the exported model and actual
solver definitions after reloading. It does not reconstruct that input's full
historical annotation environment. The optional full annotation-pipeline flag
is `--energy-candidate E5` with `--offline --no-solve` and a new output path;
that broader pipeline still requires its separately configured research inputs.

The push was prepared from committed code plus the energy-specific changes.
Other pending GPR, gene annotation and lipid work was excluded. See
[push_verification.json](push_verification.json) for the checks on that exact
export, separately from the earlier repair results.

## File identities

- Input SHA256: `aad701126d12d113816fda4b872333b614b4469ee1b6d8ab8c419231c89e965f`
- E5 SHA256: `43042d16f874a91d61f12f9f2b65838cbc264821d98bcf6e108300bf787007be`

E0 preserves loaded model definitions, not necessarily original XML bytes.
COBRA's export makes twenty formerly absent tRNA-residue charge attributes
explicitly zero. This reflects its loading convention; it is not independent
chemical evidence. Original input bytes are preserved.
