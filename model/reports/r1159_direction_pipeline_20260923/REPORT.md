# R1159 conditional leak direction integrated

User authorization on 2026-09-23: modify the reviewed [-1000,0] direction, integrate the pipeline, and push. Stored stoichiometry remains H+[cy] → H+[go]; permitted net transport is Golgi → cytosol. This assumes an acidic Golgi and no reversal of the electrochemical driving force by membrane potential. W29/PO1f pH, membrane potential and leak mechanism remain unverified; 1000 is inherited model capacity, not measured permeability. Reviewed sources and locators are versioned in `data/reference_build/curation/r1159_direction.json`.

The default builder applies a guarded curation after final field selection. Unexpected identity, chemistry, compartment, bounds, GPR or conflicting notes raises an error before edits. Original stoichiometry, GPR and all other model structure are preserved. The existing provenance collector includes the new code and curation.

Validation: local targeted regression passes (1 test); isolated staged-source regression passes (16 tests, including reaction selection and CoQ9 checks). Two full offline/no-solve builds completed in 70.36 and 66.23 seconds. Independent verification compared the complete XML trees, allowed only R1159 upper bound and exact evidence notes, and verified all recorded input/source/output hashes. No LP/screen/FVA was run.

The local new model is `model_metadata_trna_r1159_leak.xml`, compared to the fixed R539-labelled default model: 2314 reactions, 1877 metabolites, 1073 gene entries. The isolated staged build is a different source baseline and was compared to tracked `model_metadata_trna.xml`; it is not the screened V-ATPase candidate and does not promote its test-only GPR. Prior models and screen results remain unchanged. The saved R1159-zero WT witness explains why this bound alone does not establish essentiality; there is no new essentiality screen result.

Publication scope: exactly five pipeline/curation/test/README files. The local model and detailed build artifacts are retained locally; pending unrelated source, curation and model changes are not included in this push. Complete identities and dirty-state preservation checks are in `before.json`, the two build manifests and `validation.json`.

First test attempts were interrupted while unittest generated a huge textual diff: the existing semantics helper formats solver variables differently after model.copy (0 versus 0.0, signed zero). The expected model now reads the same fixture directly; full comparison and export assertions remain. Original test/logs and concrete formatting differences are preserved in the attempt1 files. No biological threshold or constraint was relaxed.

Independent audit by r1159_integration_audit reviewed exact staged scope, integration placement, sign and preconditions, opened both actual outputs and test logs, and independently repeated whole-XML and hash comparisons. All implementation checks passed; the species-specific driving force remains a scientific limitation, not an experimentally confirmed property. Commit/push receipt is recorded separately after completion.

Push completed: commit `f21e5dc5f75dd30cafb499dcded55da10566c7f8` on `origin/codex/r989-gpr-main-worktree`; exact remote ref verified. Published content matches the isolated tested source. See `push_record.json`.
