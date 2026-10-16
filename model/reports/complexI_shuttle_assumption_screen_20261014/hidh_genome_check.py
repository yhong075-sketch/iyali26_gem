"""Which N-terminus of homoisocitrate DH does the W29 genome encode? Reads the saved NC_090773.1 region.

Region: NC_090773.1:1323290-1325646 (gene YALI1_D13238g at 1323890-1325047, 1-based; NCBI Gene 2910431 esummary gives 0-based chrstart 1323889), saved under
$IYALI26_RESEARCH_ROOT/artifacts/tasks/complexI_shuttle_assumption_screen_20261014/hidh/region.fa.
Run: python3 -I model/reports/complexI_shuttle_assumption_screen_20261014/hidh_genome_check.py
Writes hidh_genome_check.json here. Never overwrites.
"""
import hashlib
import itertools
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGION = Path(os.environ["IYALI26_RESEARCH_ROOT"]) / "artifacts/tasks" / HERE.name / "hidh/region.fa"
REGION_SHA = "47af56fad06005a33a3ccccf28349b3455e35bd0190c0131aac9bcea17fcdea8"
BASES = "TCAG"
CODE = {"".join(c): aa for c, aa in zip(itertools.product(BASES, repeat=3),
                                        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
EXON1_END = "MFANRAMRLQTRMYSTGAKTLKI"  # exon 1 ends with the first base of the G codon
EXON2_START = "LIPGDGIGREVIPAG"


def translate(seq):
    return "".join(CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def main():
    assert hashlib.sha256(REGION.read_bytes()).hexdigest() == REGION_SHA
    s = "".join(line.strip() for line in REGION.read_text().splitlines() if not line.startswith(">")).upper()
    frame0 = translate(s)
    start = frame0.find(EXON1_END) * 3
    assert start >= 0
    donor = start + 3 * len(EXON1_END) + 1      # after the first base of the G codon
    exon2 = frame0.find("NPG" + EXON2_START) * 3 + 7  # 'GT' completes the G codon; then CTT (L)
    intron = s[donor:exon2]
    spliced = s[start:donor] + s[exon2:]
    protein = translate(spliced).split("*")[0]
    upstream = frame0[:start // 3]
    result = {
        "region_sha256": REGION_SHA,
        "donor_context": s[donor - 10:donor] + "|" + s[donor:donor + 10],
        "acceptor_context": s[exon2 - 20:exon2] + "|" + s[exon2:exon2 + 12],
        "intron_length": len(intron), "intron_5prime": intron[:6], "intron_3prime": intron[-3:],
        "branch_TACTAAC_offset": intron.find("TACTAAC"),
        "spliced_protein_length": len(protein), "spliced_protein_nterm": protein[:30],
        "upstream_peptide_after_last_stop_no_M": upstream.rsplit("*", 1)[-1],
        "model_start_MPPQ_inside_intron": donor < frame0.find("MPPQILLDSATEKW") * 3 < exon2,
    }
    fasta = (HERE / "hidh_isoforms.fasta").read_text().split(">")[1:]
    seqs = {block.splitlines()[0]: "".join(block.splitlines()[1:]) for block in fasta}
    result["spliced_matches_XP_502662.2"] = protein == seqs["HIDH_XP_502662.2"]
    result["spliced_matches_XP_502662.1_from_M13"] = protein[12:] == seqs["HIDH_XP_502662.1"]
    (HERE / "hidh_genome_check.json").open("x").write(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
