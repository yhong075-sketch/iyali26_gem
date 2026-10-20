"""Step-1 retrieval with Europe PMC and UniProt mocked through httpx.MockTransport."""

import asyncio
import json
import tempfile
import unittest
from pathlib import Path

import httpx
from cobra import Model, Reaction

from gem_review.lit.chunks import fulltext_chunks
from gem_review.lit.client import CachedClient, DiskCache, EPMC_FULLTEXT_URL, EPMC_SEARCH_URL
from gem_review.lit.pairs import list_pairs, select_reactions
from gem_review.lit.retrieve import retrieve_all, write_outputs
from gem_review.lit.terms import build_terms, clean_reaction_name, locus_spellings

JATS = """<article><front/><body>
<sec><title>Results</title>
<p>Unrelated paragraph about mannitol.</p>
<p>Deletion of <italic>POX2</italic> (YALI0F10857g) abolished activity on C18 acyl-CoA.</p>
<p>EC 1.3.3.6 activity was measured in vitro.</p>
</sec></body></article>"""


def _model():
    model = Model("t")
    rxn = Reaction("R101")
    rxn.name = "acyl-CoA oxidase (octadecanoyl-CoA)"
    rxn.subsystem = "Fatty acid degradation"
    rxn.annotation["ec-code"] = "1.3.3.6"
    rxn.gene_reaction_rule = "YALI1D32754g"
    other = Reaction("R999")
    other.subsystem = "Purine metabolism"
    other.gene_reaction_rule = "YALI1A00001g"
    nogpr = Reaction("R500")
    nogpr.subsystem = "Sphingolipid metabolism"
    model.add_reactions([rxn, other, nogpr])
    model.genes.get_by_id("YALI1D32754g").annotation["uniprot"] = "Q9P8Q9"
    return model


class MockApi:
    """Canned Europe PMC / UniProt answers; counts every network call."""

    def __init__(self):
        self.calls = []

    def handler(self, request: httpx.Request) -> httpx.Response:
        self.calls.append(str(request.url))
        url = str(request.url).split("?")[0]
        if url == EPMC_SEARCH_URL:
            query = request.url.params["query"]
            if "1.3.3.6" in query:
                return httpx.Response(200, json={"hitCount": 2, "resultList": {"result": [
                    {"pmid": "111", "pmcid": "PMC111", "isOpenAccess": "Y",
                     "abstractText": "POX2 encodes an acyl-CoA oxidase."},
                    {"pmid": "222", "pmcid": None, "isOpenAccess": "N",
                     "abstractText": "Lipid production in Yarrowia."}]}})
            return httpx.Response(200, json={"hitCount": 40, "resultList": {"result": []}})
        if url == EPMC_FULLTEXT_URL.format(pmcid="PMC111"):
            return httpx.Response(200, text=JATS)
        if url.startswith("https://rest.uniprot.org/uniprotkb/"):
            return httpx.Response(200, json={"genes": [{"geneName": {"value": "POX2"},
                                                       "synonyms": [{"value": "AOX2"}],
                                                       "orfNames": [{"value": "YALI0F10857g"}]}]})
        return httpx.Response(404)


def _run(api: MockApi, cache_dir: Path, pairs):
    client = CachedClient(DiskCache(cache_dir), transport=httpx.MockTransport(api.handler),
                          min_interval=0.0)

    async def go():
        try:
            return await retrieve_all(client, pairs, counterpart=lambda g: "yali0f10857g",
                                      uniprot_accession=lambda g: "Q9P8Q9")
        finally:
            await client.aclose()

    return asyncio.run(go()), client


class PairSelectionTest(unittest.TestCase):
    def test_refuses_whole_model(self):
        with self.assertRaises(ValueError):
            select_reactions(_model())

    def test_tier_and_explicit_ids(self):
        model = _model()
        self.assertEqual([r.id for r in select_reactions(model, tier=1)], ["R101", "R500"])
        self.assertEqual([r.id for r in select_reactions(model, reaction_ids=["R999"])], ["R999"])
        with self.assertRaises(KeyError):
            select_reactions(model, reaction_ids=["R_missing"])
        pairs = list_pairs(select_reactions(model, tier=1), 1)
        self.assertEqual([(p.gene, p.reaction_id, p.ec_codes) for p in pairs],
                         [("YALI1D32754g", "R101", ("1.3.3.6",))])


class TermsTest(unittest.TestCase):
    def test_queries(self):
        self.assertEqual(locus_spellings("YALI1D32754g"), ["YALI1D32754g", "YALI1_D32754g"])
        # The crosswalk hands back lower-case YALI0 tags.
        self.assertEqual(locus_spellings("yali0d24750g"), ["YALI0D24750g", "YALI0_D24750g"])
        self.assertEqual(clean_reaction_name("acyl-CoA oxidase (octadecanoyl-CoA)"), "acyl-CoA oxidase")
        self.assertIsNone(clean_reaction_name("R12"))
        terms = build_terms("YALI1D32754g", "YALI0F10857g", ["POX2", "YALI0F10857g"],
                            "acyl-CoA oxidase (octadecanoyl-CoA)", ("1.3.3.6",))
        self.assertEqual(terms.gene_terms, ("YALI1D32754g", "YALI1_D32754g", "YALI0F10857g",
                                            "YALI0_F10857g", "POX2"))
        self.assertEqual(terms.pair_query(),
                         '("YALI1D32754g" OR "YALI1_D32754g" OR "YALI0F10857g" OR "YALI0_F10857g" '
                         'OR "POX2") AND "Yarrowia lipolytica" AND ("1.3.3.6" OR "acyl-CoA oxidase")')
        self.assertIsNone(build_terms("YALI1D32754g", None, [], "", ()).pair_query())


class ChunkTest(unittest.TestCase):
    def test_only_matching_paragraphs(self):
        chunks = fulltext_chunks(JATS, pmid="111", pmcid="PMC111", terms=["POX2", "1.3.3.6"])
        self.assertEqual([(c.chunk_id, c.matched_terms) for c in chunks],
                         [("PMC111:p1", ("POX2",)), ("PMC111:p2", ("1.3.3.6",))])
        self.assertIn("abolished activity", chunks[0].text)


class RetrievalTest(unittest.TestCase):
    def test_end_to_end_and_cache(self):
        model = _model()
        pairs = list_pairs(select_reactions(model, tier=1), 1)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            api = MockApi()
            results, client = _run(api, tmp / "cache", pairs)
            self.assertEqual(len(results), 1)
            r = results[0]
            self.assertEqual((r.gene_hits, r.pair_hits, r.pair_oa_hits), (40, 2, 1))
            self.assertEqual(r.symbols, ["POX2", "AOX2"])
            self.assertEqual(r.yali0, "YALI0F10857g")
            self.assertEqual(r.pmids, ["111", "222"])
            self.assertEqual(r.fulltext_fetched, ["PMC111"])
            self.assertEqual([c.chunk_id for c in r.chunks],
                             ["111:abstract", "222:abstract", "PMC111:p1", "PMC111:p2"])
            self.assertEqual(client.network_requests, 4)  # uniprot, gene, pair, fulltext

            table = write_outputs(results, out_dir=tmp / "out", report_dir=tmp / "report",
                                  manifest={"model_sha256": "x"})
            rows = table.read_text().splitlines()
            self.assertEqual(len(rows), 2)
            self.assertIn("R101,acyl-CoA oxidase (octadecanoyl-CoA),1.3.3.6,", rows[1])
            self.assertIn(",40,2,1,1,Y,2,2,111;222", rows[1])
            stored = json.loads((tmp / "out" / "YALI1D32754g__R101.json").read_text())
            self.assertTrue(stored["fulltext_available"])

            # Second run: everything comes from the disk cache, no network call.
            api2 = MockApi()
            results2, client2 = _run(api2, tmp / "cache", pairs)
            self.assertEqual(api2.calls, [])
            self.assertEqual(client2.cache_hits, 4)
            self.assertEqual(results2[0].to_dict(), r.to_dict())


if __name__ == "__main__":
    unittest.main()
