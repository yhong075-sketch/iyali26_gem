import pathlib, re, cobra
p=pathlib.Path("/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-essential-gene-prediction-e76be3/0a4ec4c9-d926-4768-b75c-c14c6fe4ecec/scratchpad/vat/m059.xml")
m=cobra.io.read_sbml_model(p)
for met in m.metabolites:
    n=(met.name or "")
    if re.search(r"\biron\b|Fe2\+|Fe3\+|^Fe|copper|calcium|zinc|^cysteine|L-cysteine$",n,re.I) or (met.formula or "") in ("Fe","Cu","Ca","Zn"):
        print(met.id,repr(n),met.compartment,met.formula,met.charge,[(r.id,r.bounds) for r in met.reactions][:6])
