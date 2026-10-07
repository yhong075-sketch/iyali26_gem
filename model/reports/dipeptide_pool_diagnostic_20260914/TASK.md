# Four-dipeptide pool: bounded vacuolar-network diagnostic

Authorization: user's 2026-09-14 request “如果先加上一个四种二肽的pool 测试一下液泡网络能不能修复”. This authorizes a temporary supply hypothesis and targeted calculation, not acceptance as native chemistry or a permanent model/GPR/medium change.

Question: does restoring access to the four dipeptides permit vacuolar hydrolysis, water transport and R795 proton pumping? Does any restored pumping become required for maximum growth?

Input: reuse the exact three-common-AND hypothesis model, SD-Leu medium, PO1f overlay and recorded loader code from artifacts/r795_r1363_open_diagnostic_20260911/results.json after full SHA checks. Later R1025/R1026 development models are not substituted. Historical whole dirty environment is not reconstructed; record current loaded dependencies and protected files.

Pool interpretation: four independent, irreversible artificial source reactions into the existing CYTOSOL species: Gly-Asp m1870[C_cy], Gly-Glu m1861[C_cy], Ala-Gly m1877[C_cy], Gly-Pro m1865[C_cy]. Each has bounds [0,1] in model flux units. These are diagnostic boundary supplies without a GPR or a claimed enzyme; they do not interconvert the four chemicals and do not assert endogenous synthesis. The capacity is a declared test value, not a measured concentration, finite initial stock, or physiological production rate. All original exchange uptake bounds stay unchanged. The existing cytosol-to-vacuole transports remain part of the test.

Open R795 and R1363 to [0,1000] in both baseline and pool contexts, following the earlier diagnostic. No other pre-existing bounds, chemistry, GPRs or objectives change except explicitly listed per-solve objectives, R795-off controls and the near-optimal growth floor.

Exactly 12 primary LPs, 1 CPU, 60 seconds each, 600 seconds outer limit, one attempt with no automatic retry or parameter adjustment:
1. No pool, maximum growth.
2. No pool, maximum R795, original biomass lower bound.
3. Pool, maximum growth.
4. Pool, maximum R795, original biomass lower bound.
5. Pool, maximum R1363, original biomass lower bound.
6–9. Pool, separately maximize R2021/R2029/R2034/R2039 hydrolysis.
10. Pool, R795 off, maximum growth.
11. Pool, minimum R795 with biomass >= pool WT minus 1e-8 h^-1.
12. Pool, R795 off, maximize the sum of the four hydrolyses.

All runs retain original directionality and all other connections. This is a reaction-level control, not a gene knockout screen or a new essentiality-count benchmark. Explicit predictions outside these tests remain unknown. Existing chemistry/formula limitations are not silently repaired.

Acceptance: optimal finite results; full steady-state and bound residual <=1e-7; baseline growth matches the prior value within abs 1e-8 / rel 1e-6. Positive flux means >1e-7 for this diagnostic, not a gene-essentiality threshold. Save all fluxes and bounds, exact artificial sources, complete relevant compartment species rows, full run input/source identities and actual solver settings. Independently audit equations, directions, computed results and the distinction between restored feasibility and biological repair. Stop on identity conflict, failure/nonfinite output, or exhausted budget; preserve partial evidence. No model export, production curation change, full screen, cluster work, git commit/push, or physical experiment.
