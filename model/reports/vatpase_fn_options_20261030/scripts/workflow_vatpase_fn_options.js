export const meta = {
  name: 'vatpase-fn-options',
  description: 'Read-only analysis of why V-ATPase genes are FN in iYali26 and which fixes are evidence-backed',
  phases: [
    { title: 'Gather', detail: 'history, model structure, literature, experimental calls in parallel' },
    { title: 'Design', detail: 'one agent builds the option set' },
    { title: 'Verify', detail: 'three skeptics with distinct lenses' },
  ],
}

const REPO = '/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/essential-gene-prediction-e76be3'
const RES = '/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research'
const SP = '/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-essential-gene-prediction-e76be3/0a4ec4c9-d926-4768-b75c-c14c6fe4ecec/scratchpad'
const MODEL = SP + '/vat/m059.xml'
const COMMON = `
Context: iYali26 is a cobrapy genome-scale model of Yarrowia lipolytica. In the current candidate (SHA-256 059316ef34417c88aadba61cfd5606134d27879a842b3835d2cc7df917a9248d, a copy at ${MODEL}) eleven V-ATPase genes on reactions R794/R795 are false negatives: experimentally essential (consensus list ${RES}/state/essentiality/repository/consensus_essential_genes.csv) but predicted non-essential (KO/WT about 1.0). The 15 % KO/WT cutoff is primary; Task A medium + PO1f profile. The user asked: "how can we solve the VTPase problem?"
Repo (read-only): ${REPO}. Reports in model/reports/, decision log model/STATE.md (newest first; older entries are in Chinese), evaluation rules model/expected/benchmark.md.
HARD RULES: read-only. Do not edit any repo file or model. Do NOT run any LP/FBA/FVA/optimization (no model.optimize, no slim_optimize, no solver). Static reading, cobra loading and Boolean GPR evaluation are allowed (python: /opt/miniconda3/bin/python; pass the model path as pathlib.Path). Write scratch files only under ${SP}/vat/agents/. Tag every factual claim verified (you checked it this run and show the evidence), inferred, or recalled. For literature claims, "verified" means you opened that exact URL this run and read the passage; a search summary is not a source. Never print secrets.`

const HIST_SCHEMA = { type: 'object', properties: {
  timeline: { type: 'array', items: { type: 'object', properties: { date: {type:'string'}, report: {type:'string'}, what_was_done: {type:'string'}, result: {type:'string'} }, required: ['date','report','what_was_done','result'] } },
  user_decisions: { type: 'array', items: { type: 'string' } },
  established_mechanism: { type: 'string' },
  open_questions: { type: 'array', items: { type: 'string' } },
  options_already_tried_or_rejected: { type: 'array', items: { type: 'string' } },
}, required: ['timeline','user_decisions','established_mechanism','open_questions','options_already_tried_or_rejected'] }

const MODEL_SCHEMA = { type: 'object', properties: {
  r794: { type: 'string' }, r795: { type: 'string' },
  vatpase_genes: { type: 'array', items: { type: 'object', properties: { gene: {type:'string'}, annotation: {type:'string'}, in_branch: {type:'string'}, experimental_essential: {type:'string'}, ko_closes_r794: {type:'boolean'}, ko_closes_r795: {type:'boolean'} }, required: ['gene','annotation','in_branch','experimental_essential','ko_closes_r794','ko_closes_r795'] } },
  proton_network: { type: 'string', description: 'every reaction producing/consuming/transporting H+ in vacuole and Golgi, with bounds and whether H+-gradient-coupled' },
  gradient_dependent_processes: { type: 'array', items: { type: 'string' } },
  other_bypasses: { type: 'array', items: { type: 'string' } },
  checks_run: { type: 'array', items: { type: 'string' } },
}, required: ['r794','r795','vatpase_genes','proton_network','gradient_dependent_processes','other_bypasses','checks_run'] }

const LIT_SCHEMA = { type: 'object', properties: {
  claims: { type: 'array', items: { type: 'object', properties: { claim: {type:'string'}, url: {type:'string'}, tier: {type:'string', enum:['verified','recalled','inferred','established_uncited']}, organism: {type:'string'}, quote_or_passage: {type:'string'} }, required: ['claim','url','tier','organism','quote_or_passage'] } },
  ylipolytica_direct_evidence: { type: 'string' },
  essentiality_mechanisms: { type: 'array', items: { type: 'string' } },
  how_other_gems_handle_it: { type: 'string' },
  unread_or_failed_urls: { type: 'array', items: { type: 'string' } },
}, required: ['claims','ylipolytica_direct_evidence','essentiality_mechanisms','how_other_gems_handle_it','unread_or_failed_urls'] }

const DATA_SCHEMA = { type: 'object', properties: {
  table: { type: 'array', items: { type: 'object', properties: { gene: {type:'string'}, subunit_guess: {type:'string'}, in_model: {type:'boolean'}, consensus_essential: {type:'boolean'}, raw_calls: {type:'string'} }, required: ['gene','subunit_guess','in_model','consensus_essential','raw_calls'] } },
  consistency: { type: 'string' },
  sources: { type: 'array', items: { type: 'string' } },
}, required: ['table','consistency','sources'] }

phase('Gather')
const [hist, model, lit, data] = await parallel([
  () => agent(`${COMMON}
TASK: Reconstruct the full history of the V-ATPase non-essentiality problem in this repo. Read every report folder whose name contains vatpase, r794_r795, r795_open, r795_r1363, vacuole_open_supply, iyli647_vatpase, screen_vatpase, um2_vacuolar_closures_20261028, plus model/reports/essentiality_fn_causes_20261027/REPORT.md, and every model/STATE.md entry mentioning V-ATPase, R794, R795, R1159 or UM2 (grep for them). Report what was tried, each result with numbers, the user's recorded decisions (quote them), the established mechanism, and open questions. Also list options already tried or rejected and why.`, { label: 'history', phase: 'Gather', schema: HIST_SCHEMA }),
  () => agent(`${COMMON}
TASK: Static and Boolean analysis of the current model ${MODEL}.
1. Print R794 and R795: equation with compartments, bounds, full GPR, notes.
2. For each gene in either GPR: annotation from the model, which branch, whether it is in the consensus essential list, and whether a single knockout makes each GPR false (cobra GPR eval).
3. List every reaction that produces, consumes or transports H+ in the vacuole and Golgi compartments (find the compartment IDs first), with bounds and whether each couples H+ movement to another solute (antiport/symport) or ATP. Use the stoichiometric matrix, not names.
4. From that, list the processes in the model that would need a vacuolar/Golgi proton gradient (H+-coupled transporters) and whether each is open, closed, or blocked; note whether any of them is required for biomass (static reasoning only: is the transported metabolite in the biomass reaction, and are there other routes?).
5. List any other route by which vacuolar/Golgi H+ balance can be satisfied without R794 (e.g. R1159 diffusion), with bounds.
Also read ${SP}/next/current_fn_fp.json and report which V-ATPase-related genes are listed as FN and with what cause label.`, { label: 'model-structure', phase: 'Gather', schema: MODEL_SCHEMA }),
  () => agent(`${COMMON}
TASK: Literature. Answer with sources you actually open (use WebFetch / WebSearch; prefer Europe PMC, PubMed, journal full text):
(a) Is there DIRECT experimental evidence that V-ATPase subunits (VMA genes, VPH1/STV1) are essential or required in Yarrowia lipolytica? Include genome-wide CRISPR screens (e.g. Schwartz et al. 2019, Patterson et al. 2018 Tn-seq, Robertson 2024) and any targeted mutant.
(b) Why are V-ATPase mutants lethal or conditionally lethal in fungi? In S. cerevisiae vma mutants are viable at acidic pH but fail at pH 7.5, on non-fermentable carbon, with high Ca2+/metals; find primary sources. Find the mechanism for respiratory failure (e.g. vacuolar amino acid / cysteine storage and iron homeostasis, Hughes et al.), and evidence for essentiality in Crabtree-negative or obligately aerobic fungi (Candida albicans, Kluyveromyces lactis, Neurospora crassa, Aspergillus, Magnaporthe).
(c) How do genome-scale models (Yeast-GEM / Yeast8, iYali, human Recon, others) represent V-ATPase, and has any published method captured V-ATPase essentiality (pH constraints, membrane potential, ME-models, maintenance terms)?
Give one claim per row with URL, tier and the passage read. List URLs you could not open.`, { label: 'literature', phase: 'Gather', schema: LIT_SCHEMA }),
  () => agent(`${COMMON}
TASK: Experimental data consistency. Identify every Y. lipolytica gene that is a V-ATPase subunit or assembly factor (V1 A-H, V0 a/c/c'/c''/d/e, VMA21, VMA12/VPH2, VMA22, PKR1, VOA1, RAV1/RAV2) using the model's gene annotations, the repo's earlier V-ATPase reports (model/reports/*vatpase*, r794_r795_gpr_review_20260909) and the consensus list ${RES}/state/essentiality/repository/consensus_essential_genes.csv. Also look for the raw screen calls under ${RES}/state/essentiality/ (any per-screen files with Cas9/Cas12a/Tn-seq scores). For each gene give: likely subunit, in model or not, consensus essential or not, and raw calls with scores. Then say whether the data are consistent (all subunits essential?) and whether anything looks like a mapping or call problem. Keep raw values; do not swap labels.`, { label: 'experimental-calls', phase: 'Gather', schema: DATA_SCHEMA }),
])

phase('Design')
const OPT_SCHEMA = { type: 'object', properties: {
  diagnosis: { type: 'string' },
  options: { type: 'array', items: { type: 'object', properties: {
    id: {type:'string'}, name: {type:'string'}, what_changes: {type:'string'}, mechanism_and_evidence: {type:'string'},
    expected_effect_on_fn: {type:'string'}, conflicts_with_decisions_or_rules: {type:'string'}, evidence_still_missing: {type:'string'},
    authorization_and_compute_needed: {type:'string'}, recommendation: {type:'string'} },
    required: ['id','name','what_changes','mechanism_and_evidence','expected_effect_on_fn','conflicts_with_decisions_or_rules','evidence_still_missing','authorization_and_compute_needed','recommendation'] } },
  recommended_path: { type: 'string' },
}, required: ['diagnosis','options','recommended_path'] }

const gathered = JSON.stringify({ history: hist, model: model, literature: lit, experimental: data }, null, 1)
const design = await agent(`${COMMON}
You have four read-only findings (JSON below). Build the option set for solving the V-ATPase false negatives.
Rules from the project's AGENTS.md that bind every option: never change the model, a GPR, a bound or biomass merely to improve recall; adding a biomass demand or closing a bypass needs mechanistic evidence, ideally direct Y. lipolytica evidence; respect recorded user decisions (e.g. UM2 Option A0 of 2026-10-28) and say explicitly when an option would reopen one; model changes go through curated data + build; computation needs explicit authorization with a budget.
Cover at least: (1) GPR correction alone; (2) GPR + an evidence-based gradient-coupled requirement (name which process, from the model findings and literature, and whether the model can represent it); (3) a generic acidification/maintenance pseudo-demand; (4) reclassifying these genes as outside what a stoichiometric model can predict (a reporting class, no model change); (5) anything else the findings suggest. For each option state the expected effect on these 11 FN as inferred unless a check proves it, and what would make it fail (e.g. the second OR branch, R1159, other H+ routes). Then give one recommended path.
FINDINGS:
${gathered}`, { label: 'design', phase: 'Design', schema: OPT_SCHEMA })

phase('Verify')
const VERDICT = { type: 'object', properties: {
  per_option: { type: 'array', items: { type: 'object', properties: { id: {type:'string'}, verdict: {type:'string', enum:['holds','partly_holds','refuted','unverifiable']}, reason: {type:'string'}, evidence: {type:'string'} }, required: ['id','verdict','reason','evidence'] } },
  errors_in_findings: { type: 'array', items: { type: 'string' } },
  missing: { type: 'array', items: { type: 'string' } },
}, required: ['per_option','errors_in_findings','missing'] }
const designJson = JSON.stringify(design, null, 1)
const LENSES = [
  ['rules', 'Lens: project rules and recorded decisions. Read AGENTS.md and model/STATE.md yourself. Try to show that an option is recall tuning, reopens a recorded user decision without saying so, needs computation it does not name, or touches the colleague-owned lipid-unlump line. Also check the history agent quotes the decisions correctly.'],
  ['mechanics', `Lens: model mechanics. Load ${MODEL} yourself (no LP). Try to refute each option's expected effect: would the 11 genes really become essential, or does the second OR branch, R1159, another H+ route, a closed reaction, or a blocked transporter defeat it? Re-run the Boolean GPR checks and stoichiometric traces that the claims depend on and show the output.`],
  ['evidence', 'Lens: biological evidence. Open every URL the literature findings tag verified and check that the passage supports the exact claim and organism. Downgrade anything unread. Try to refute the claim that the V-ATPase is essential in Y. lipolytica, and the claimed mechanism linking V-ATPase loss to growth failure.'],
]
const verdicts = await parallel(LENSES.map(([k, lens]) => () => agent(`${COMMON}
You are an adversarial reviewer. ${lens}
Default to refuted or unverifiable when you cannot check a claim yourself.
OPTION SET:
${designJson}
UNDERLYING FINDINGS:
${gathered}`, { label: 'skeptic:' + k, phase: 'Verify', schema: VERDICT })))

return { hist, model, lit, data, design, verdicts: { rules: verdicts[0], mechanics: verdicts[1], evidence: verdicts[2] } }