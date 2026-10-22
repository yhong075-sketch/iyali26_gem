"""Prompt caching: where the breakpoint lands and the order calls run in (Anthropic API mocked).

The request bodies are captured from the real pydantic-ai -> anthropic SDK path through an httpx2
mock transport, so these tests check the bytes that would be sent, not our own description of them.
"""

import asyncio
import json
import os
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("PYDANTIC_AI_NO_BANNER", "1")

import httpx2  # noqa: E402
from anthropic import AsyncAnthropic  # noqa: E402
from pydantic_ai.models.anthropic import AnthropicModel  # noqa: E402
from pydantic_ai.providers.anthropic import AnthropicProvider  # noqa: E402

from gem_review.lit.extract import build_agent, run_extraction  # noqa: E402
from gem_review.lit.providers import settings_for  # noqa: E402

ANTHROPIC = "anthropic:claude-opus-5-5"
# Long enough that the chunk block, not the rules text, is clearly the bulk of the prefix.
TEXT_A = "Aox3 was found to be active against short-chain fatty acids in Y. lipolytica. " * 40
TEXT_B = "Dga1p acylates diacylglycerol in the endoplasmic reticulum of Y. lipolytica. " * 40


def _record(gene, rid, name, chunk_text, pmid):
    chunk = {"chunk_id": f"{pmid}:abstract", "pmid": pmid, "source": "abstract", "text": chunk_text}
    return {"gene": gene, "reaction_id": rid, "reaction_name": name, "ec_codes": ["1.3.3.6"],
            "yali0": None, "symbols": [], "subsystem": "x", "tier": 1, "pair_hits": 3,
            "pmids": [pmid], "fulltext_available": False, "chunks": [chunk]}


SHARED = [_record("YALI1D32754g", f"R{n}", f"acyl-CoA oxidase (C{n})", TEXT_A, "111") for n in (4, 6, 8)]
ONE_OFF = _record("YALI1E38810g", "R1771", "diacylglycerol acyltransferase", TEXT_B, "222")


class Recorder:
    """Mock Messages API: records each request body and the order requests start and finish."""

    def __init__(self):
        self.bodies: list[dict] = []
        self.events: list[tuple[str, str]] = []

    async def handler(self, request):
        body = json.loads(request.content)
        question = body["messages"][0]["content"]
        question = question if isinstance(question, str) else question[0]["text"]
        rid = question.split("reaction ")[1].split(" ")[0]
        self.bodies.append(body)
        self.events.append(("start", rid))
        await asyncio.sleep(0.01)
        self.events.append(("end", rid))
        return httpx2.Response(200, json={
            "id": "msg", "type": "message", "role": "assistant", "model": "claude-opus-5-5",
            "content": [{"type": "text", "text": '{"claims": []}'}], "stop_reason": "end_turn",
            "stop_sequence": None,
            "usage": {"input_tokens": 50, "output_tokens": 5, "cache_creation_input_tokens": 900,
                      "cache_read_input_tokens": 0}})


def _run(records, model_name=ANTHROPIC):
    rec = Recorder()
    client = AsyncAnthropic(api_key="test", http_client=httpx2.AsyncClient(transport=httpx2.MockTransport(rec.handler)))
    model = AnthropicModel("claude-opus-5-5", provider=AnthropicProvider(anthropic_client=client),
                           settings=settings_for(ANTHROPIC))
    with tempfile.TemporaryDirectory() as tmp:
        results = asyncio.run(run_extraction(records, agent=build_agent(model), model_name=model_name,
                                             out_dir=Path(tmp)))
        usage = [json.loads(line) for line in (Path(tmp) / "usage.jsonl").read_text().splitlines()]
    return rec, results, usage


def _by_reaction(rec):
    out = {}
    for body in rec.bodies:
        question = body["messages"][0]["content"]
        question = question if isinstance(question, str) else question[0]["text"]
        out[question.split("reaction ")[1].split(" ")[0]] = body
    return out


class BreakpointTest(unittest.TestCase):
    def test_breakpoint_is_on_the_chunk_block_not_the_rules_or_question(self):
        rec, _, _ = _run(SHARED)
        for body in rec.bodies:
            rules, chunks = body["system"]
            self.assertTrue(rules["text"].startswith("You review the evidence"))
            self.assertNotIn("cache_control", rules)
            self.assertTrue(chunks["text"].startswith("Retrieved chunks (JSON):"))
            self.assertEqual(chunks["cache_control"], {"type": "ephemeral", "ttl": "5m"})
            user = body["messages"][0]["content"]
            blocks = [{"text": user}] if isinstance(user, str) else user
            self.assertFalse(any("cache_control" in b for b in blocks))   # the question stays uncached
            self.assertNotIn("cache_control", body)                       # no automatic tail caching
            self.assertEqual(body["output_config"]["effort"], "medium")   # agent settings kept
            self.assertEqual(body["max_tokens"], 16000)

    def test_prefix_is_byte_identical_across_reactions_sharing_chunks(self):
        rec, _, _ = _run(SHARED)
        prefixes = {json.dumps({k: v for k, v in b.items() if k != "messages"}, sort_keys=True)
                    for b in rec.bodies}
        self.assertEqual(len(rec.bodies), 3)
        self.assertEqual(len(prefixes), 1)

    def test_one_off_chunk_set_is_not_cached(self):
        rec, results, usage = _run([ONE_OFF, *SHARED])
        body = _by_reaction(rec)["R1771"]
        self.assertFalse(any("cache_control" in block for block in body["system"]))
        self.assertEqual({u["reaction_id"]: u["cached"] for u in usage},
                         {"R1771": False, "R4": True, "R6": True, "R8": True})

    def test_other_providers_get_no_cache_settings(self):
        rec, _, usage = _run(SHARED, model_name="openai:gpt-5")
        self.assertFalse(any("cache_control" in block for body in rec.bodies for block in body["system"]))


class OrderAndUsageTest(unittest.TestCase):
    def test_first_call_of_a_group_finishes_before_the_others_start(self):
        rec, _, _ = _run([ONE_OFF, *SHARED])
        shared = [(kind, rid) for kind, rid in rec.events if rid in {"R4", "R6", "R8"}]
        first = shared[0][1]
        self.assertEqual(shared[1], ("end", first))   # nothing else in the group starts before it ends
        self.assertEqual(sorted(rid for kind, rid in shared[2:] if kind == "start"),
                         sorted({"R4", "R6", "R8"} - {first}))

    def test_cache_tokens_are_logged(self):
        _, results, usage = _run(SHARED)
        self.assertTrue(all(r.status == "ok" for r in results))
        for row in usage:
            self.assertEqual((row["cache_write_tokens"], row["cache_read_tokens"]), (900, 0))
            self.assertEqual(row["input_tokens"], 50 + 900)  # pydantic-ai counts cached tokens in input_tokens


if __name__ == "__main__":
    unittest.main()
