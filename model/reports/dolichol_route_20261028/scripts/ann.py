import cobra, warnings, sys
warnings.simplefilter("ignore")
M = cobra.io.read_sbml_model(sys.argv[1])
for rid in ("R2177","R2178","R2179","R2180","R2181","R2182","R287","R1101","R288","R2106","R2105","R2138","R217","R235","R230","R231"):
    r = M.reactions.get_by_id(rid)
    print(rid, r.name, "| notes:", dict(r.notes), "| ann:", {k:v for k,v in r.annotation.items() if k!='sbo'}, "| subsystem:", r.subsystem)
for mid in ("m2019[C_cy]","m2018[C_cy]","m499[C_cy]","m501[C_cy]","m503[C_er]","m412[C_lp]","m411[C_lp]","m432[C_lp]","m1969[C_er]","m1968[C_er]"):
    m = M.metabolites.get_by_id(mid)
    print(mid, m.name, dict(m.notes), {k:v for k,v in m.annotation.items() if k!='sbo' and k!='inchi'})
