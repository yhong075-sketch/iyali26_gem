# Evidence files (fetched or extracted 2026-10-29)

| File | Source | Note |
| --- | --- | --- |
| `uniprot_aars_atp_synthase_20261029.tsv` | UniProt REST, `rest.uniprot.org/uniprotkb/search?query=accession:…&fields=accession,id,gene_orf,gene_oln,gene_names,xref_geneid,protein_name,length,ft_transit,cc_subcellular_location,sequence&format=tsv` for the 16 accessions named in REPORT.md §3 | SHA-256 `303f78b11266494879d14cb377304a746fff4ef52d03688f0bf75e1448faa305`; fetched with curl; CLIB122 (`YALI0`) entries mapped to the model's `YALI1` genes through `$IYALI26_RESEARCH_ROOT/reference/locus_map/s2_metabolic_genes.csv` |
| `cas9_cas12a_calls_target_genes.csv` | Rows of `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/cas9_cas12a_fitness.csv` (SHA-256 `ba1eca8f…` source workbook; normalized table) for the 22 genes of the five rules | Genes missing from the file were not in either screen's table (YALI1D21476g and the three mtDNA genes) |
