# Independent figure evidence audit

Reviewed at 2026-09-21T23:15:44.277618+00:00 by the independent `branch_function_review` subagent. **The displayed reaction data, saved results and scientific wording pass this bounded review. The revised WT and shared-component KO PNG exports have also been directly inspected; the initial label collisions were corrected before acceptance.** No optimizer, new prediction or scientific model mutation was used by this audit.

## Direct checks

The audit directly parsed the screened candidate XML, stored runtime bounds and full WT/target-knockout flux files. It separately read the generated `map.json`, `evidence.json`, HTML template and provenance. It did not treat generator assertions as independent verification.

- All nine source files listed by the figure provenance retain their full SHA-256 identities.
- R794, R795 and R1159 retain every reactant/product coefficient in the XML. All three displayed GPRs match the XML Boolean AST exactly; R1159 has no assigned GPR. Runtime WT bounds are R794=[0,1000], R795=[0,0], R1159=[−1000,1000]. Displayed reversibility matches these bounds.
- Every one of the 15 saved scenarios (WT plus 14 target single KOs) matches its stored growth, KO/WT ratio, disabled GPRs, newly closed reactions and all three displayed fluxes. All 14 systematic IDs are preserved.
- For the 12 shared candidate knockouts, R794 is newly closed and R795 was already closed. The two a-candidate single KOs keep both GPRs true. The earlier independent screen audit verified the complete WT flux remains feasible under every target KO; this figure introduces no new optimization.
- The original WT has zero flux through all three drawn reactions. R1159 transfers H⁺ between cytosol and Golgi reversibly with no ATP term in its stored stoichiometry. Its presence supports a structural route for review, not a claim that an active bypass explains the current WT.
- `biomass_C` consumes ATP but has no direct Golgi-H⁺ or vacuolar-H⁺ term. This is a statement about that reaction and the saved screen, not a claim that every indirect compartment connection is absent.

The local adjacency read found five reactions involving Golgi H⁺ and 22 involving vacuolar H⁺. All have zero saved WT flux. Of the latter, six have bounds permitting nonzero flux and 16 are closed, including R795 and R1163. The figure intentionally draws only three core reactions and explicitly says that other proton-linked reactions and the wider ATP network are omitted. Nonzero bounds do not establish network feasibility or actual activity.

## Claim ledger

“Supported” below means the limited displayed computational statement matches the inspected evidence. “Unresolved” identifies the scientific question explicitly left open by the figure; it is not a failed data check or evidence that the native function is absent.

| # | Displayed statement or open question | Review status |
|---|---|---|
| 1 | The three drawn reactions preserve current candidate stoichiometry, GPR and bounds. | Supported |
| 2 | The 15 displayed scenarios reproduce saved results rather than new live optimization. | Supported |
| 3 | Shared-component KO propagation works, while the saved WT flux remains feasible and growth stays approximately unchanged. | Supported; current model and condition only |
| 4 | R1159 is a reversible proton connection with no ATP term and zero saved WT flux. | Supported; potential structural alternative, not observed active bypass |
| 5 | Biomass ATP demand alone does not impose required pump flux in this screen. | Supported by the saved feasible WT and direct biomass stoichiometry |
| 6 | Which measured, condition-specific acidification requirement should constrain growth? | Unresolved |
| 7 | What physiological direction and energetic coupling correctly describe R1159? | Unresolved |
| 8 | Can the two native a candidates support the same compartment under the same conditions? | Unresolved |
| 9 | Which intact, expressed W29 protein supplies the old F candidate's function? | Unresolved |
| 10 | What is the native complete assembly and which functional dependencies are confirmed? | Unresolved |

**10 total / 10 audited / 5 supported / 5 unresolved / 0 contradicted / 0 unchecked.** The figure's five evidence-gap headings correctly retain the unresolved questions; their existence must not be rephrased as five proven biological defects.

## Biological and structural scope

Gene functions remain candidates and native established names are unverified. The figure preserves the prior review's distinctions: YALI1E12482g (V0 a candidate, native localization unresolved); YALI1F38820g (Vph1-like V0 a candidate, experimental confirmation required); YALI0E16192g (V1 F candidate represented by a CLIB122 sequence, W29 identity unresolved). The complete gene/function table is present in the interactive HTML. These identities come from the existing source/sequence/GPR review, which this audit reused rather than conducting a new primary-literature investigation.

The one-copy AlphaFold branch inputs do not establish a native stoichiometric assembly, catalysis, localization or redundancy. Previous branch A raw files were audited; branch B remained webpage-summary-only in the reviewed evidence. The HTML states this asymmetry and does not claim the prediction comparison resolved native function.

The wording does not prescribe an arbitrary acidification sink, ATP demand increase, direction change or GPR strengthening to force an essential call. It distinguishes the verified computational propagation from the unverified native biology. All plotted reactions and their detailed panels are partial-model views, not a newly accepted model.

## Artifact identity

- Candidate XML SHA: `0fcc2f6ff05124b91977c6f96f78cbe5906a3ccbf3f936f95990afe02b540353`.
- Map SHA at data audit: `3c79398ff57f8d978df1a1c22ee6c00e8b1ea494fdb6d5fd6dd37825f7dea224`.
- Evidence JSON SHA at data audit: `9be34b53efd73c42ce925df10a266de79cc52cfea550670bc7ead525cfbead75`.
- Full upstream source identities are in `provenance.json`; screen-level independent checks are in `../screen_vatpase_candidate_20260921/AUDIT.md` and its `audit_checks.json`.

## Visual and final-data review

The auditor directly viewed the actual WT and shared-component KO PNG exports. The initial WT image had two label collisions and an ambiguous spacing-only a-candidate caption. After the layout revision, the H₂O/ATP/ADP labels are unobscured; the two a-candidate descriptions are separate; the four coefficient-2 proton labels are visible; the KO image shows R794 as newly closed (red dashed) and R795 as already closed (gray dashed). The selected KO is clearly identified as YALI1D00581g, a candidate V1 D central-shaft subunit with unverified native name/function. The figure preserves “all drawn WT fluxes = 0,” stored-directionality, partial-core-view and test-only-GPR qualifications. No material scientific overclaim or remaining visual obstruction was found in the reviewed revision.

A small presentation recommendation was sent to include the primary KO/WT <15% threshold beside the static essential-gene count. The root agent's `render-check.json` reports all 15 browser scenarios, drag/edit and export checks passed with no page errors; this is a root-executed UI check, not an independently repeated browser interaction. This auditor independently repeated the data binding checks after the layout revision (3 reactions, 15 saved scenarios and 9 protected source hashes), without running a solver.

## Entry-page repair: independent static check (2026-09-21T23:27:27.254036+00:00)

The root agent reported that the user's actual Chrome tab opened `viewer.html`. This auditor directly read that file before repair: it was the build template, with `const data=__DATA__` still present, so normal JavaScript evaluation would stop before populating the selector or constructing the map. This supports the reported blank-page symptom; the actual tab inventory was the root agent's observation, not an independently repeated browser inspection.

The fix reaches that same filename. `viewer.html` now contains a zero-second relative meta refresh plus a clickable fallback link to the existing `vatpase-escher-en.html`. Both resolve to an existing file. The original template is preserved byte-for-byte as `viewer.html.in`; the builder reads that filename and checks that neither directly openable `.html` file contains the placeholder.

An independent before/after SHA comparison confirmed byte-identical `map.json`, `evidence.json`, the generated `vatpase-escher-en.html`, `escher.min.js`, and all four WT/KO SVG/PNG exports. All nine upstream source hashes still match. Only the entry/template naming and build plumbing changed; the prior scientific and visual audit conclusions are unchanged. This audit performed static inspection and hashing only. No new browser reproduction was performed or claimed; the root agent reported that CUA access to the local-file tab was blocked by its URL safety policy and stopped that access attempt.
