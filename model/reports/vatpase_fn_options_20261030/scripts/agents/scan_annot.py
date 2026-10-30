import openpyxl, re, sys
R="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"
wb=openpyxl.load_workbook(R+"/reference/locus_map/S2_table_YALI1_YALI0_map.xlsx",read_only=True)
pat=re.compile(r"\b(VMA\d+|VPH1|VPH2|STV1|CUP5|TFP1|VAT2|TFP3|PPA1|PKR1|VOA1|RAV1|RAV2|YDL185W|YBR127C|YKL080W|YEL051W|YOR332W|YGR020C|YHR039C-A|YHR039C-B|YPR036W|YOR270C|YMR054W|YEL027W|YPL234C|YHR026W|YLR447C|YCL005W-A|YGR105W|YKL119C|YHR060W|YMR123W|YGR106C|YJR033C|YDR202C)\b|vacuolar.{0,20}ATP|H\+-ATPase|V-ATPase|V0|V1 ", re.I)
for ws in wb.worksheets:
    for i,row in enumerate(ws.iter_rows(values_only=True),start=1):
        txt=" ".join(str(x) for x in row if x is not None)
        if pat.search(txt):
            print(ws.title,i,row)
