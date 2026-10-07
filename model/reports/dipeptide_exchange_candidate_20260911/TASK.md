# Four-exchange diagnostic supply candidate

User request: fix the Exchange section first. Clarification offered: test-model uptake versus diagram-only. Proceeding with the stated test-model interpretation; preserve original scientific inputs.

Persistent change: a separate medium CSV, SD-Leu plus Gly-Asp/Gly-Glu/Ala-Gly/Gly-Pro, each uptake cap1000 (diagnostic only, not measured). Load it through the existing effective-simulation-context loader so media application cannot silently reset these exchanges. Expected effective exchange bounds [-1000,1000]; preserve all other medium entries and the PO1f profile. Do not edit production/default medium or original XML/GPR/stoichiometry.

Context: same common-AND hypothesis XML; R795 and R1363 temporarily [0,1000] in both old/new comparison contexts. No exported scientific model. Controls: original medium versus supplemented diagnostic medium.

Run exactly8 primaryLPs: old WT, candidate WT, candidate maximum R795, candidate maximum R1363, and candidate minimum of each of4 exchange reactions. Flux objectives retain original biomass lower bound0, without demanding optimal growth. No KO or full screen. Save all2315 fluxes and actual bounds per solve, raw finite status/values, identity/version/runtime parameters and supply assumptions.

Budget600seconds outer, Threads1, TimeLimit60/LP, presolveFalse, no retry. Stop on mismatched identities, nonoptimal/nonfinite results or numerical balance/bound violation>1e-7. Old/current genotype overlay operation effects and non-target reaction bounds must match. Independent source audit checks actual input/output, not claims alone. Display any updated Escher as a new diagnostic scenario with its own provenance; retain the prior figure and results.

## Current status after the user asked what the dipeptides are

Prepared only; not executed. No optimization, knockout, model export, default-medium change or figure update has occurred. The separate CSV and runner are a proposed exogenous-supply diagnostic, not an established correction to the experimental medium. Whether to activate external dipeptide supply or investigate internal production remains unresolved. Existing original sources and original medium hashes are unchanged. Do not describe this candidate as an applied Exchange fix.
