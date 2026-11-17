"""Step 2, static part (no LP): the matrix CoA-moiety balance of every reaction.

Moiety = the pantetheine unit: species whose name contains "CoA", "coenzyme A" or "pantetheine" (CoA, its
thioesters, 3'-dephospho-CoA, pantetheine 4'-phosphate, pantetheine). The species list is written out so it can be
checked by hand. For each reaction, delta_mi = sum of the stoichiometric coefficients of moiety species in C_mi:
- 0 and every moiety species in C_mi: internal, conserving;
- 0 with moiety species inside and outside C_mi: a pass-through transport (moves the moiety, net 0 into C_mi);
- not 0: changes the matrix moiety pool (a boundary transport, or a matrix synthesis or consumption step).
Writes ../tables/matrix_moiety_reactions.csv, ../tables/moiety_species.csv and
../tables/coa_species_reactions.csv (producers, consumers, transporters per species of interest).
"""
import csv
import re
import warnings
from pathlib import Path

import cobra

from common import MODEL, check_inputs

OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(exist_ok=True)
PATTERN = re.compile(r"CoA|coenzyme A|pantetheine", re.I)
# Species asked about in the task: CoA and every acyl-CoA that can cross between C_mi and C_cy.
FOCUS = ["m73[C_mi]", "m69[C_cy]", "m78[C_mi]", "m1875[C_cy]", "m1970[C_cy]", "m156[C_mi]", "m2046[C_cy]",
         "m2047[C_cy]", "m648[C_mi]", "m646[C_cy]", "m1831[C_mi]", "m2050[C_cy]", "m1828[C_mi]", "m2049[C_cy]",
         "m1829[C_mi]", "m1971[C_cy]", "m71[C_mi]", "m68[C_cy]", "m475[C_mi]", "m1350[C_mi]"]

check_inputs()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = cobra.io.read_sbml_model(MODEL)
moiety = {m.id for m in model.metabolites if PATTERN.search(m.name or "")}

with open(OUT / "moiety_species.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["species", "name", "formula", "charge", "compartment", "n_reactions"])
    for m in sorted(model.metabolites, key=lambda x: (x.compartment, x.id)):
        if m.id in moiety and m.compartment in ("C_mi", "C_cy"):
            w.writerow([m.id, m.name, m.formula, m.charge, m.compartment, len(m.reactions)])


def kind(r):
    mi = {m: c for m, c in r.metabolites.items() if m.id in moiety and m.compartment == "C_mi"}
    out = {m: c for m, c in r.metabolites.items() if m.id in moiety and m.compartment != "C_mi"}
    delta = sum(mi.values())
    if not mi:
        return None, None
    if delta == 0 and not out:
        return delta, "internal (conserving)"
    if delta == 0:
        return delta, "pass-through (net 0 into C_mi)"
    if out:
        return delta, "boundary transport"
    return delta, "matrix source or sink"


rows = []
for r in model.reactions:
    delta, k = kind(r)
    if k is None:
        continue
    rows.append([r.id, r.name, r.lower_bound, r.upper_bound, r.gene_reaction_rule or "", delta, k,
                 r.build_reaction_string()])
rows.sort(key=lambda x: (x[6], x[0]))
with open(OUT / "matrix_moiety_reactions.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["reaction", "name", "lb", "ub", "gene_rule", "delta_moiety_into_C_mi_per_unit_forward_flux",
                "class", "equation"])
    w.writerows(rows)
counts = {}
for row in rows:
    counts[row[6]] = counts.get(row[6], 0) + 1
print("reactions touching the matrix moiety:", len(rows), counts)
for row in rows:
    if row[6] != "internal (conserving)":
        print(row[6], "|", row[0], row[1], (row[2], row[3]), repr(row[4]), "delta", row[5], "|", row[7])

with open(OUT / "coa_species_reactions.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["species", "name", "reaction", "role", "lb", "ub", "gene_rule", "other_compartments", "equation"])
    for sid in FOCUS:
        m = model.metabolites.get_by_id(sid)
        for r in sorted(m.reactions, key=lambda x: x.id):
            c = r.metabolites[m]
            comps = sorted({x.compartment for x in r.metabolites} - {m.compartment})
            if comps and any(x.compartment != m.compartment and x.id in moiety for x in r.metabolites):
                role = "transport"
            elif r.lower_bound < 0 < r.upper_bound:
                role = "producer and consumer (reversible)"
            elif (c > 0 and r.upper_bound > 0) or (c < 0 and r.lower_bound < 0):
                role = "producer"
            elif r.lower_bound == 0 and r.upper_bound == 0:
                role = "closed"
            else:
                role = "consumer"
            w.writerow([sid, m.name, r.id, role, r.lower_bound, r.upper_bound, r.gene_reaction_rule or "",
                        ";".join(comps), r.build_reaction_string()])
print("wrote", OUT)
