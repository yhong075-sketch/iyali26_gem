"""In-memory test of proposed acyl-ACP formula fills (nothing saved). Prints every touched reaction's residual."""
import json, sys, cobra
M = cobra.io.read_sbml_model(sys.argv[1])
FILLS = {  # acyl group + S (pantetheine thiol) + R (rest of ACP), as for acetyl-ACP C2H3OSR and octanoyl-ACP C8H15OSR
 "m1410[C_mi]": "C4H5O2SR",   # acetoacetyl
 "m1415[C_mi]": "C4H7O2SR",   # 3-hydroxybutanoyl
 "m1418[C_mi]": "C4H5OSR",    # trans-but-2-enoyl
 "m1411[C_mi]": "C4H7OSR",    # butanoyl
 "m1412[C_mi]": "C6H9O2SR",   # 3-oxo-hexanoyl
 "m1416[C_mi]": "C6H11O2SR",  # 3-hydroxyhexanoyl
 "m1419[C_mi]": "C6H9OSR",    # trans-hex-2-enoyl
 "m1413[C_mi]": "C6H11OSR",   # hexanoyl
 "m1414[C_mi]": "C8H13O2SR",  # 3-oxo-octanoyl
 "m1417[C_mi]": "C8H15O2SR",  # 3-hydroxyoctanoyl
 "m1420[C_mi]": "C8H13OSR",   # trans-oct-2-enoyl
 "m1901[C_mi]": "C10H19OSR", "m1910[C_cy]": "C10H19OSR",   # decanoyl
 "m1902[C_mi]": "C12H23OSR", "m1911[C_cy]": "C12H23OSR",   # dodecanoyl
 "m1903[C_mi]": "C14H27OSR", "m1912[C_cy]": "C14H27OSR",   # myristoyl
 "m1904[C_mi]": "C14H25OSR", "m1913[C_cy]": "C14H25OSR",   # cis-tetradec-7-enoyl
 "m1905[C_mi]": "C16H31OSR", "m1914[C_cy]": "C16H31OSR",   # palmitoyl
 "m1906[C_mi]": "C16H29OSR", "m1915[C_cy]": "C16H29OSR",   # cis-hexadec-9-enoyl
 "m1907[C_mi]": "C18H35OSR", "m1916[C_cy]": "C18H35OSR",   # octadecanoyl
 "m1908[C_mi]": "C18H33OSR", "m1917[C_cy]": "C18H33OSR",   # cis-octadec-11-enoyl
 "m1909[C_mi]": "C18H31OSR", "m1918[C_cy]": "C18H31OSR",   # "octadecynoyl" (C18H31O acyl; product named octadecadienoate)
 "m1919[C_cy]": "C8H15OSR",   # octanoyl (cytosol)
 "m1920[C_cy]": "HSR",        # cytosolic ACP
 "m1922[C_cy]": "C18H32O2",   # octadecadienoate
}
touched = set()
for mid, f in FILLS.items():
    m = M.metabolites.get_by_id(mid); print(f"{mid:13s} {m.name[:40]:40s} {str(m.formula):10s} -> {f}"); m.formula = f
    touched |= {r.id for r in m.reactions}
print()
for rid in sorted(touched, key=lambda x: int(x[1:])):
    r = M.reactions.get_by_id(rid)
    b = {k: round(v, 3) for k, v in r.check_mass_balance().items() if abs(v) > 1e-9}
    nadph = sum(c for m, c in r.metabolites.items() if m.name.startswith("NADPH"))
    print(f"{rid:6s} {str(r.bounds):18s} {b}  (NADPH coef {nadph:g}) | {r.build_reaction_string(use_metabolite_names=False)[:110]}")
