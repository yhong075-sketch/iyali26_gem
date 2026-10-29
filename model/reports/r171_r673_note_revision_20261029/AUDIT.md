# Independent audit: R171 / R673 note revision, 25dabf62 vs 6af03c40 (2026-10-29)

A fresh read-only general-purpose subagent with no conversation history audited this round. It edited nothing and ran
no LPs. It re-read Hahn 2016 (Europe PMC PMC4980432) and the live UniProt REST records. Its scratch scripts were in
the session scratchpad and are not kept. Summary of its report, with the values it returned:

1. **SHAs and build records: supported.**
   - Models: old `6af03c40…`, new `25dabf62…`.
   - `before/` copies: `63a7c6ba…` (R171 data), `2446217a…` (R673 data), `315646ab…` (step).
   - Current files: `1b1d719f…`, `e1d151a6…`, `cf395b25…`.
   - The seed-0 and seed-1 `.build.json` files differ only in the output path, timestamps and dirty status; both give
     `25dabf62…`.
   - Against the old candidate's build record, only the two data SHAs and the step source SHA changed; HEAD
     `ff7fc62`.
   - Re-running both `write_curation.py` scripts in a scratch tree reproduces the curated files byte for byte.
2. **XML diff: supported.**
   - Raw diff: 2 changed `<p>` lines (70620, 94953), identical to `logs/raw_xml_diff.txt`.
   - cobra semantic diff: only R171 `r171_atp_synthase_gpr_evidence` and R673 `r673_prolyl_split_evidence`.
   - Both new values match TASK §3 character for character.
3. **Curated JSONs vs `before/`: supported.**
   - R171: only the evidence note changed.
   - R673: the evidence note, the P38708 and P39965 claims and the new P78600 entry. These are exactly the fields
     TASK §3 lists.
   - The docstring changed by one line.
4. **Wording vs evidence: supported, with one minor imprecision and one candidate overstatement.**
   - R171 is fixed. Hahn 2016 says the monomer "lacks subunits e, g, and k", and gives ATP hydrolysis of both
     monomer and dimer at about 2.25 U/mg. The R171 audit's other gaps (1) and (2) stay in REPORT §1 / §9, not in
     the note.
   - R673 is fixed. P38708 is "Putative…", with no subcellular-location comment; P39965 is "Probable…,
     mitochondrial", ECO:0000269, PMID 16823961; P78600 is "…cytoplasmic", ECO:0000305. The annotation strings
     confirm both clauses.
   - **Imprecision (verified):** P38708 carries the GO CC term cytoplasm (IEA:InterPro). So "no UniProt location"
     holds only for the Subcellular location section. The curated `evidence[2]` wording, "no subcellular location
     annotated", is the precise one.
   - **Candidate overstatement (inferred):** "the assignment rests on the two predictors and synteny". The
     "ancestral locus Anc_x" tags come from the same "similar to *S. cerevisiae* …" annotation clause. As far as the
     auditor knows, *Y. lipolytica* is not in YGOB (established / uncited). So this "synteny" is not independent of
     the homology names, and no *Y. lipolytica* gene-order comparison is shown.
5. **Optimisation problem unchanged: supported (structurally verified).**
   - The stoichiometric matrices, the reaction / metabolite / gene order, the bounds, the GPRs and the objective are
     identical.
   - The solve path (`essentiality_simulation_context.py`) does not read notes.
   - Sanity solves as recorded: WT 0.3116608216; `YALI1F04433g` 0.0341900 (0.10970, R171 closed); `YALI1E05929g` 0
     (R673 closed).

## Defects and actions

| # | Defect | Action |
| --- | --- | --- |
| 1 | Full pytest still running | Done after the audit (REPORT §5) |
| 2 | No `build_commands.txt`; seed not in `.build.json` | `logs/build_commands.txt` added; the seed gap is inherent to the builder (disclosed) |
| 3 | `notes_and_solves.json` lacks medium, profile and solver provenance | Script extended and re-run (3 more LPs, same values); it now records the medium `6eca70fd…`, the profile `35307853…` and Gurobi at 1e-9 |
| 4 | "no UniProt location" imprecise; "predictors and synteny" may overstate | **Not changed**: it would change the wording the user approved. Recorded in REPORT §6 and put to the user |
| 5 | libSBML `R_R304` notes error in the build logs | Pre-existing (every recent build); not this round's |
| 6 | REPORT.md missing at audit time | Written |

## Second pass (v3 `059316ef` vs `25dabf62`)

A fresh read-only subagent audited the second pass; no edits, no LPs; scratch not kept.
1. **SHAs and build records: supported.** v3 `059316ef…`, seed 1 identical; data `cdb9dbc0…` (before `e1d151a6…`);
   against the v2 record only the R673 data SHA changed; all `source_sha256` equal; `input_unchanged` and
   `source_unchanged_during_build` true.
2. **XML diff: supported.** One hunk, 4 lines; ElementTree and cobra diffs: only R673 notes
   `r673_prolyl_split_evidence` and `r673_prolyl_split_evidence_tier`; applying the TASK §6 replacements to the old
   values reproduces the new ones exactly.
3. **Data and regeneration: supported.** Exactly the three listed fields changed; the generator in a scratch copy
   reproduces `cdb9dbc0…` byte for byte.
4. **Wording: supported** for both changes. "No UniProt subcellular-location annotation" is precise (P38708 has an
   empty Subcellular location section; GO CC cytoplasm IEA and multienzyme complex IBA); the ancestral-locus tags sit
   in the same "similar to S. cerevisiae …" clause. Residual (from pass 1): P78600 "annotated cytoplasmic" holds for
   its name and Subcellular location section, but its GO CC is mitochondrion (IEA:TreeGrafter; re-fetched by the main
   session, [evidence/uniprot_P78600_go_cc.tsv](evidence/uniprot_P78600_go_cc.tsv)).
5. **Remaining in the R673 data: partially.** `evidence[3]` "synteny partner" (disclosed); `evidence[2]` "best
   S. cerevisiae hit" (no search shown; the source is the annotation string); model note
   `r673_prolyl_split_previous_gpr` and `unlinked_gene.identity` say "AIM10-like", one of the two mixed names;
   `evidence[*].tier` "ortholog reference" overstates "similar to".

| # | Defect | Action |
| --- | --- | --- |
| 1–4 | Residual wording (best hit; synteny partner; P78600 GO; AIM10-like; ortholog tier) | **Not changed** (outside the approved fields); recorded in REPORT §8 as limits |
| 5 | No pass-2 build_commands; no copy of the pre-pass-2 generator; GO fields not saved | `pass2/logs/build_commands.txt` added; GO CC fields saved (`evidence/*_go_cc.tsv`); the generator's pass-2 edit is exactly the TASK §6 replacements (no file copy kept) |
| 6 | Full pytest running | See REPORT §8 |
