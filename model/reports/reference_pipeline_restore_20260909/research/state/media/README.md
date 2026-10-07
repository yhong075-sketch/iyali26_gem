# Experimental media

`sd_leu.csv` models the screen condition reported by Ramesh et al. (2023):
synthetic defined medium lacking leucine, with 2% glucose, at 30 °C.

At runtime the validator closes uptake through every exchange reaction, then
opens only the exchanges listed in this file. `R1219` (L-leucine exchange) is
intentionally absent and must remain closed.

The finite amino-acid/nucleobase bounds are concentration-scaled estimates:
glucose is 111 mM at uptake 10 mmol/gDW/h, and each CSM-Leu component is scaled
relative to that reference. Concentration is not an uptake rate, so these values
are assumptions rather than fitted parameters. The validator reports cutoff
sensitivity without fitting these bounds to the essential-gene labels.

The 0.67% YNB fraction also supplies 0.2 mg/L FeCl3.6H2O, equivalent to
0.000740 mM Fe. Because iYali26 exposes only the `R1189` iron(2+) exchange,
the medium uses that reaction as a lumped bioavailable-iron boundary and gives
it a concentration-scaled uptake bound of 0.0000667 mmol/gDW/h. The Fe(III) to
Fe(II) speciation step is therefore outside the model boundary. This bound
records experimental availability; it is not fitted to essentiality recall.

Sources:

- https://doi.org/10.1038/s42003-023-04996-8
- https://tools.thermofisher.com/content/sfs/manuals/IFU459942.pdf
