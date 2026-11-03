"""W29 tag-protein audit: flag logic, coverage of every model gene, decisions for every disagreement, and a rebuild
from the raw research-workspace inputs; no LPs, no cobra."""

import json
import unittest

from scripts.gem_annotate.locus_resolver import model_gene_fingerprint
from tools.audit_w29_tag_protein import (
    FLAGS,
    NEEDS_DECISION,
    SIDECAR_NAME,
    audit_from_research_root,
    audit_path,
    build_audit,
    check_decisions,
    classify,
    model_gene_ids,
    read_audit,
    read_decisions,
)


def s2_row(strand="+", chrom="YALI1A", start=100, end=400, gene_type="mRNA"):
    return dict(yali0="YALI0A00001g", gene_type=gene_type, chrom=chrom, start=start, end=end, strand=strand)


def cds_row(strand="+", chrom="YALI1A", start=100, end=400):
    return dict(protein="XP_000001.1", chrom=chrom, start=start, end=end, strand=strand)


class ClassifyTests(unittest.TestCase):
    def test_flags(self):
        gene = "YALI1A00001g"
        cases = {
            "ok": (gene, s2_row(), cds_row(start=350, end=900)),
            "strand_flip": (gene, s2_row(), cds_row(strand="-")),
            "no_overlap": (gene, s2_row(), cds_row(start=401, end=900)),
            "no_s2_row": (gene, None, cds_row()),
            "no_refseq_cds": (gene, s2_row(), None),
            "mito_absent_from_feature_table": ("YALI1M00001g", s2_row(chrom="YALI1M"), None),
            "not_protein_coding": ("YALI1A00001r", s2_row(gene_type="tRNA"), None),
            "not_yali1": ("YALI0A00001g", None, None),
        }
        self.assertEqual(set(cases), set(FLAGS))
        for flag, args in cases.items():
            with self.subTest(flag=flag):
                self.assertEqual(classify(*args), flag)

    def test_other_chromosome_is_no_overlap_and_strand_wins(self):
        self.assertEqual(classify("YALI1A00001g", s2_row(), cds_row(chrom="YALI1B")), "no_overlap")
        self.assertEqual(classify("YALI1A00001g", s2_row(), cds_row(strand="-", start=5000, end=6000)), "strand_flip")

    def test_non_yali1_ids(self):
        for gene in ("YALIUNK1", "YalifMp05", "YALIfMp29", "YALI0E16192g"):
            with self.subTest(gene=gene):
                self.assertEqual(classify(gene, s2_row(), cds_row()), "not_yali1")

    def test_non_mrna_row_with_a_cds_is_still_compared(self):
        self.assertEqual(classify("YALI1A00001t", s2_row(gene_type="RetroElement"), cds_row(strand="-")),
                         "strand_flip")


class CommittedAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = read_audit()
        cls.sidecar = json.loads(audit_path(SIDECAR_NAME).read_text())
        cls.decisions = read_decisions()

    def test_audit_covers_exactly_the_model_genes(self):
        genes = model_gene_ids()
        self.assertEqual({row["gene"] for row in self.rows}, genes,
                         "model genes changed: rerun python -m tools.audit_w29_tag_protein")
        self.assertEqual(len(self.rows), len(genes))
        self.assertEqual(self.sidecar["model_gene_count"], len(genes))
        self.assertEqual(self.sidecar["model_gene_fingerprint"], model_gene_fingerprint(genes))

    def test_flags_are_known_and_match_the_sidecar(self):
        counts = {flag: sum(row["flag"] == flag for row in self.rows) for flag in FLAGS}
        self.assertEqual({row["flag"] for row in self.rows} - set(FLAGS), set())
        self.assertEqual(counts, self.sidecar["flag_counts"])

    def test_every_disagreement_has_a_decision(self):
        self.assertEqual(check_decisions(self.rows, self.decisions), [])
        self.assertTrue(any(row["flag"] in NEEDS_DECISION for row in self.rows))


GENE = "YALI1A00001g"


def flipped_audit():
    return build_audit([GENE], {"YALI1_A00001g": s2_row()}, {"YALI1_A00001g": cds_row(strand="-")})


def decision(**changes):
    row = dict(schema_version="1", gene=GENE, flag="strand_flip", refseq_protein="XP_000001.1",
               refseq_span="YALI1A:100-400-", decision="refseq_protein_is_antisense_orf", sequence_source="UniProt",
               sequence_accession="Q00001", evidence="run 1", reason="antisense")
    row.update(changes)
    return row


class CheckDecisionsTests(unittest.TestCase):
    def test_valid_decision_passes(self):
        self.assertEqual(check_decisions(flipped_audit(), [decision()]), [])
        self.assertEqual(check_decisions(flipped_audit(), [decision(decision="unresolved_conflict",
                                                                    sequence_source="none",
                                                                    sequence_accession="")]), [])

    def test_each_rule_reports(self):
        cases = {
            "has no decision": [],
            "duplicate decision": [decision(), decision()],
            "unsupported schema_version": [decision(schema_version="2")],
            "is not allowed for flag": [decision(decision="refseq_protein_kept", sequence_source="RefSeq",
                                                 sequence_accession="XP_000001.1")],
            "unknown sequence_source": [decision(sequence_source="GenBank")],
            "sequence_accession is required": [decision(sequence_accession="")],
            "must point to the CLIB122 UniProt protein": [decision(sequence_source="RefSeq",
                                                                   sequence_accession="XP_000001.1")],
            "must not name a protein": [decision(decision="unresolved_conflict")],
            "evidence and reason are required": [decision(reason="")],
            "decision is for flag": [decision(flag="no_overlap", decision="refseq_protein_kept",
                                              sequence_source="RefSeq", sequence_accession="XP_000001.1")],
            "decision was made for": [decision(refseq_protein="XP_000002.1")],
        }
        for message, decisions in cases.items():
            with self.subTest(message=message):
                self.assertTrue(any(message in error for error in check_decisions(flipped_audit(), decisions)),
                                check_decisions(flipped_audit(), decisions))

    def test_kept_refseq_protein_must_be_the_audited_one(self):
        audit = build_audit([GENE], {"YALI1_A00001g": s2_row()}, {"YALI1_A00001g": cds_row(start=500, end=900)})
        kept = decision(flag="no_overlap", refseq_span="YALI1A:500-900+", decision="refseq_protein_kept",
                        sequence_source="RefSeq", sequence_accession="XP_000001.1")
        self.assertEqual(check_decisions(audit, [kept]), [])
        errors = check_decisions(audit, [dict(kept, sequence_accession="XP_000009.1")])
        self.assertIn(f"{GENE}: a kept RefSeq protein must name that RefSeq protein", errors)

    def test_resolved_disagreement_leaves_a_stale_decision(self):
        audit = [dict(flipped_audit()[0], flag="ok")]
        self.assertIn(f"{GENE}: stale decision, the audit no longer flags this gene",
                      check_decisions(audit, [decision()]))


class RebuildTests(unittest.TestCase):
    def test_raw_inputs_reproduce_the_committed_audit(self):
        # Fails closed: load_project_paths(required=True) raises when the research workspace is not configured.
        rows, sidecar = audit_from_research_root()
        self.assertEqual(rows, read_audit())
        self.assertEqual(sidecar, json.loads(audit_path(SIDECAR_NAME).read_text()))


if __name__ == "__main__":
    unittest.main()
