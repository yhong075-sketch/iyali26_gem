# R673 prolyl-tRNA synthetase pair: TargetP 2.0 and DeepLoc 2.1 predictions (2026-10-29)

Scope and authorization: [TASK.md](TASK.md) (user, 2026-10-29: "Yes, run predictors"). Prediction stage only: no
model, GPR, curated data, code or test changed; nothing committed. The proposed rule change is TASK.md §4, waiting for
approval.

## 1. Result

Both predictors agree with the synteny assignment of the OR-rule review (§7.2): `YALI1E05929g` cytosolic,
`YALI1D19110g` mitochondrial. These are **predictions**, not localisation evidence; no *Y. lipolytica* experiment was
found for either protein.

| Gene (W29 sequence) | TargetP 2.0 (Non-Plant) | DeepLoc 2.1 (High-quality, ProtT5) |
| --- | --- | --- |
| `YALI1E05929g` (XP_503566.2, 528 aa) | **Other** 0.999947; SP 0.000020; mTP 0.000033; no cleavage site | **Cytoplasm** 0.7524 (threshold 0.4761); Mitochondrion 0.2215 (0.6220); Nucleus 0.2950 (0.5014); soluble 0.7707; signal label: nuclear localization signal |
| `YALI1D19110g` (XP_502869.1, 445 aa) | **mTP** 0.625112; Other 0.374754; cleavage 13–14 (`RRL-KD`, Pr 0.3944) | **Mitochondrion** 0.8443 (0.6220); Cytoplasm 0.2491; Nucleus 0.2081; soluble 0.7760; signal label: mitochondrial transit peptide |

Table: [localization_predictions.tsv](localization_predictions.tsv). Raw outputs: [raw/targetp/](raw/targetp/) (job
`6AE3BC76002509BB05EB4B1C`) and [raw/deeploc/](raw/deeploc/) (job `6AE3BCA7002509D23808033B`); DeepLoc's two attention
PNGs (0.75 and 0.96 MB) are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r673_prolyl_localization_20261029/deeploc/`,
SHA-256 in `raw/deeploc/attention_png_sha256.txt`. DeepLoc thresholds were read from the job's `results.json`.

The decision rule fixed in TASK.md §3 before the run is met: `YALI1E05929g` has no mitochondrial call from either
tool, and `YALI1D19110g` has one from both.

## 2. Verification (verified this session)

[scripts/verify_predictions.py](scripts/verify_predictions.py) → [prediction_verification.json](prediction_verification.json),
all checks passed:
- both jobs report 0 failed sequences and exactly the two submitted names;
- DeepLoc's attention CSVs spell each full submitted sequence (528 and 445 residues);
- TargetP's per-residue files spell the first 200 residues of each sequence (the file covers only 200, so it does not
  prove the full input on its own), equal the files in the job zip, and the zip CRCs pass.

Submission: one job per tool, both sequences together, through the DTU web forms in the built-in browser; no
e-mail address, no account. Form settings checked before submitting (TargetP organism `non-pl`, format `long`;
DeepLoc encode `Slow`, format `long`). The DeepLoc page reported ProtT5 on a GPU.

## 3. What the predictions mean and do not mean

- **Annotation context (added after the curation-stage audit).** The OR-rule review's synteny reading is not the
  whole annotation: `YALI1D19110g`'s annotation string also reads "similar to uniprot|P78600 *Candida albicans*
  Prolyl-tRNA synthetase cytoplasmic", and YHR020W (`YALI1E05929g`'s hit) is "Putative proline--tRNA ligase" with no
  UniProt location. The homology names alone are therefore mixed; the predictions above are the main
  compartment evidence. The "synteny assignment" in §1 is the annotation's ancestral-locus tag, part of the same
  homology clause, so it is not independent of those names (inferred; added by the note revision). Details: [../r673_prolyl_split_20261029/REPORT.md](../r673_prolyl_split_20261029/REPORT.md) §1.

- TargetP and DeepLoc share sequence information; the two agreeing calls are not two independent experiments, and
  their scores are not probabilities to multiply.
- `YALI1D19110g`'s TargetP mTP likelihood is moderate (0.625; cleavage site Pr 0.39); DeepLoc's mitochondrial score
  is well above its threshold (0.844 vs 0.622). Both rest on the same N-terminus `MLPRLSRLTTRRLKDV…`.
- `YALI1E05929g`: "Other" in TargetP means no recognised N-terminal signal, not "cytosol confirmed". DeepLoc's NLS
  label and nucleus score 0.295 (below 0.5014) are noted; a nuclear pool is not excluded.
- **Start site.** The W29 gene model of `YALI1E05929g` starts 11 residues upstream of the CLIB122 model (Q6C6Z6,
  517 aa; TASK §2). The predictions used the W29 sequence. The 517-aa variant was not submitted; whether it would get a
  different call is untested.
- Dual targeting of one gene product to both compartments (known for some yeast aminoacyl-tRNA synthetases,
  recalled, uncited) is not detectable by these tools.

## 4. Next gate

The rule change is proposed in TASK.md §4 and needs the user's approval.
