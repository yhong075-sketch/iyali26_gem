# Lipid candidate and audit modules

These modules are optional, explicit analyses of a frozen source model. They
never activate production changes. The canonical builder and its model.xml stay
independent of the candidates.

## Aggregate TAG acyl-moiety ledger

Run `python -m scripts.plan_lipid_moiety_ledger model.xml data/lipid_combo_curation.csv --output ledger.json`.
The compiler returns a fresh, in-memory model plus a canonical JSON manifest.
Its 95 internal reactions conserve seven acyl-chain identities, 35 TAG outputs,
and 81 terminal routes. This is an aggregate chain histogram; DAG pairing,
sn position, cross-compartment transport, quantitative TAG yield, and a complete
lipid network are outside its scope. Candidate GPR coverage does not establish
experimental enzyme assignments. Composition probabilities are metadata unless
an explicit soft-prior probe is requested.

The biochemical-pH anionic tuple convention is checked against the input.
`source_coa_acyl_formula_charge_verified` is the chemistry gate. When the source
uses another convention, `chemistry_source_not_normalized` remains a blocker;
only the separate in-memory feasibility exercise uses curated tuples. The source
is never silently repaired. Every report retains `activation_ready: false`.
A successful isolated ledger probe is not an integrated regression test of the
whole GEM.

Input SHA-256 values bind the model, CSV and specification. The CLI refuses
symlink and hard-link output aliases of its inputs and writes JSON atomically.
Determinism is checked across two fresh structural builds; solver probe results
are reused for that structural comparison.

## Strict sn-position candidate

Run `python -m scripts.lp_sn12_candidate model.xml --report strict-sn.json --candidate-sbml strict-sn.xml`.
The source must match the curated SHA-256 and full model fingerprint. The
candidate has 1,134 explicit lipid states and reactions, retains template bounds
and GPRs, rewrites biomass weights without changing their totals, and closes
specified generic bypasses. Candidate SBML uses the shared deterministic,
atomic writer. The original model remains unchanged. Reusing an already marked
candidate as source is rejected.

R989 and R39 candidate assignments remain provisional. R39's mapping is
explicitly AlphaFold-prediction-only and unverified. R1521 only receives the
recorded neutral-convention formula completion. CoA biochemical-pH migration,
connected-component chemistry and cardiolipin expansion remain blocked. No
report authorizes a production activation.

## CoA and R1521 handoffs

`scripts.fatty_acyl_coa_handoff` and `scripts.r1521_current_snapshot_handoff`
accept an explicit source model and emit read-only audit reports. They bind the
source and evidence contracts by SHA-256 and fail closed on drift. CoA
normalization is available through its guarded API but is absent from the
canonical build and raises for the checked-in blocked curation.

## Frozen provenance

The accepted source for strict-sn and handoff tests is the main-branch
`model.xml`, SHA-256
`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`.
Tests can use `IYALI26_SOURCE_MODEL` to select another copy of that exact source.
A rebuilt or modified source requires a new audit, not a relaxed expected hash.

`data/iyli21.xml` is the frozen legacy raw input, SHA-256
`6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631`,
retained only because the ER VLCFA evidence contract hashes those exact bytes.
It is distinct from the canonical `data/iyali26.xml` and does not replace the
canonical starting model. The optional ER identity correction validates this
contract, changes only its recorded targets, preserves the excluded peroxisomal
copy, and rolls back on a failed postcondition. It is not silently inserted
into the canonical builder.

The two 2026-08-18 handoff reports are historical regression fixtures, not fresh
validation claims. The supplementary Rhea/MetaNetX evidence records are retained
because the handoff loaders validate their content digests.
