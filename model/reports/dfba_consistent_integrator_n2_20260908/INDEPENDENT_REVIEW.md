# Independent numerical review

Status: **conditionally accepted as an isolated toy prototype; not accepted for production or real-GEM claims**.

The reviewer verified by static inspection that the helper has no solver dependency or production caller, uses one biomass exposure for biomass and every inventory, returns the earliest depletion event, and contains exactly three deterministic examples. The requested numerical safeguards were incorporated: a mixed endpoint-event tolerance, rejection of material pre-clamp inventory deficits, explicit overflow failure with subdivision guidance, fixed-volume/nonnegative-growth/withdrawal-only scope, and a split-step invariant.

For fixed specific rates,

\[
E=B_0\frac{\operatorname{expm1}(\mu h)}{\mu},\qquad
B_1=B_0+\mu E,\qquad C_{i,1}=C_{i,0}-q_iE,
\]

with \(E=B_0h\) when \(\mu=0\), is exact. But in FBA both \(q_i\) and \(\mu\) are optimization results, so

\[
q_iB_0\frac{e^{\mu h}-1}{\mu}\le C_i
\]

is nonlinear and coupled. It must not be presented as an ordinary pre-solve LP bound. A real implementation would require an explicitly approved event-and-resolve strategy or a separately validated self-consistent solve.

The final self-check was rerun after all requested safeguards and passed. No GEM solve was performed.
