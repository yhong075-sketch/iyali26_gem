"""Write the glycerol and oleate variants of the Task A medium for the step-5 diagnostic (new files only).

Every row, value and comment comes from model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv
(Task A, SHA-256 6eca70fd...). Only two rows change per file: D-glucose R1070 is closed (uptake 0)
and the carbon source is added at the measured W29 uptake (growth_parameters.csv, rounded to
4 decimals as in the task prompt). The CSM-Leu supplement caps keep their glucose-2.414-scaled
values, as the user chose the PO1f + CSM-Leu medium for these runs.

Usage (repository root): python -B model/reports/r132_mito_adh_20261012/build_media_step5.py
"""
import csv
import hashlib
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
MEDIA = REPO / "model/conditions/media"
SOURCE = MEDIA / "po1f_csm_leu_w29uptake_20261012.csv"
SOURCE_SHA = "6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45"
DATA = ("Worland AM et al. 2024 Metab Eng 85:1-13 doi:10.1016/j.ymben.2024.06.010 (link to the data "
        "repository inferred, pending user confirmation); data: github.com/garrettroell/yarrrowia_eflux2 "
        "commit 8361e5536aa97d83d95b77cacf0f982a154f9df0 results/growth_parameters/growth_parameters.csv")
CARBON = {
    "glycerol": ("R1141", "6.2666", "glycerol; measured_W29_uptake_Worland2024: 6.2666 mmol/gDW/h "
                 "(6.26659 +/- 0.95963) in strain W29 on a defined glycerol medium; this medium is CSM-Leu (medium differs)"),
    "oleate": ("R1440", "0.7630", "oleate; measured_W29_uptake_Worland2024: 0.7630 mmol/gDW/h "
               "(0.76302 +/- 0.07524) in strain W29 on a defined oleic acid medium; this medium is CSM-Leu (medium differs)"),
}


def main():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA, "Task A medium changed"
    rows = list(csv.DictReader(SOURCE.open()))
    fields = list(rows[0])
    for name, (rid, uptake, comment) in CARBON.items():
        assert rid not in {r["exchange"] for r in rows}
        out = MEDIA / f"po1f_csm_leu_w29{name}_20261012.csv"
        with out.open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            for row in rows:
                if row["exchange"] == "R1070":
                    row = {**row, "uptake": "0",
                           "comment": f"D-glucose; closed: {name} is the sole carbon source in this variant",
                           "uptake_basis": "closed for the carbon-source swap", "evidence_tier": "condition",
                           "source": "Task B step 5 (r132_mito_adh_20261012)"}
                writer.writerow(row)
            writer.writerow({"exchange": rid, "uptake": uptake, "comment": comment,
                             "uptake_basis": "measured_W29_uptake_Worland2024", "evidence_tier": "measured, different medium",
                             "source": DATA})
        print(out.name, hashlib.sha256(out.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
