"""Provider selection, platform/.env handling and the setup wizard (no network)."""

import os
import stat
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from gem_review.lit.providers import (DEFAULT_LLM, LLM_ENV, PROVIDERS, load_env_file, provider_for,
                                      readiness, resolve_llm, settings_for, write_env_values)
from gem_review.lit.setup_cli import run_setup


class Scripted:
    """Answers prompts in order; records every prompt and output line."""

    def __init__(self, answers: list[str]):
        self.answers = list(answers)
        self.prompts: list[str] = []
        self.lines: list[str] = []

    def read(self, prompt: str) -> str:
        self.prompts.append(prompt)
        return self.answers.pop(0)

    def out(self, line: str) -> None:
        self.lines.append(line)


class EnvFileTest(unittest.TestCase):
    def test_load_sets_missing_variables_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / ".env"
            path.write_text("# comment\nexport A_KEY='quoted value'\nB_KEY=plain\n\nnot a pair\n")
            with mock.patch.dict(os.environ, {"B_KEY": "from shell"}, clear=False):
                os.environ.pop("A_KEY", None)
                self.assertEqual(load_env_file(path), ["A_KEY"])
                self.assertEqual(os.environ["A_KEY"], "quoted value")
                self.assertEqual(os.environ["B_KEY"], "from shell")
            self.assertEqual(load_env_file(Path(tmp) / "missing"), [])

    def test_write_replaces_keeps_and_locks(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / ".env"
            path.write_text("NCBI_EMAIL=me@example.org\nOPENAI_API_KEY=old\n")
            write_env_values(path, {"OPENAI_API_KEY": "new", LLM_ENV: "openai:gpt-5"})
            self.assertEqual(path.read_text(), f"NCBI_EMAIL=me@example.org\nOPENAI_API_KEY=new\n{LLM_ENV}=openai:gpt-5\n")
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o600)


class ResolutionTest(unittest.TestCase):
    def test_flag_then_env_then_default(self):
        with mock.patch.dict(os.environ, {LLM_ENV: "openai:gpt-5"}):
            self.assertEqual(resolve_llm(None), "openai:gpt-5")
            self.assertEqual(resolve_llm("google-gla:gemini-2.5-pro"), "google-gla:gemini-2.5-pro")
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(LLM_ENV, None)
            self.assertEqual(resolve_llm(None), DEFAULT_LLM)

    def test_readiness_messages(self):
        self.assertIn("not a provider:model string", readiness("gpt-5"))
        with mock.patch.dict(os.environ, {}, clear=False):
            for var in ("OPENAI_API_KEY", "ANTHROPIC_API_KEY"):
                os.environ.pop(var, None)
            self.assertIn("No key for OpenAI", readiness("openai:gpt-5"))
            self.assertIn("OPENAI_API_KEY", readiness("openai:gpt-5"))
            os.environ["OPENAI_API_KEY"] = "x"
            self.assertIsNone(readiness("openai:gpt-5"))
        with mock.patch("gem_review.lit.providers.is_installed", return_value=False):
            self.assertIn("uv add 'pydantic-ai-slim[groq]'", readiness("groq:llama-3.3-70b-versatile"))
        self.assertIsNone(readiness("ollama:llama3.1") if provider_for("ollama:x").key_vars == () else "key needed")
        self.assertIsNone(readiness("unknown-provider:model"))  # left to pydantic-ai

    def test_settings_by_provider(self):
        self.assertEqual(settings_for("anthropic:claude-opus-5-5")["anthropic_effort"], "medium")
        self.assertEqual(settings_for("openai:gpt-5")["max_tokens"], 16000)
        self.assertNotIn("anthropic_effort", settings_for("openai:gpt-5"))


class SetupWizardTest(unittest.TestCase):
    def _run(self, answers, secret, tester, env: dict | None = None):
        script = Scripted(answers)
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, env or {}, clear=False):
            for provider in PROVIDERS:
                for var in provider.key_vars:
                    os.environ.pop(var, None)
            path = Path(tmp) / ".env"
            code = run_setup(env_path=path, read=script.read, read_secret=lambda prompt: secret,
                             tester=tester, out=script.out)
            text = path.read_text() if path.exists() else None
        return code, text, script

    def test_openai_choice_tests_and_saves(self):
        tested = []
        openai_index = next(i for i, p in enumerate(PROVIDERS, start=1) if p.prefix == "openai")
        code, text, script = self._run([str(openai_index), "", "y"], "sk-test-123", lambda m: tested.append(m))
        self.assertEqual(code, 0)
        self.assertEqual(tested, ["openai:gpt-5"])
        self.assertIn(f"{LLM_ENV}=openai:gpt-5\n", text)
        self.assertIn("OPENAI_API_KEY=sk-test-123\n", text)
        self.assertTrue(any("OPENAI_API_KEY=<hidden>" in line for line in script.lines))
        self.assertFalse(any("sk-test-123" in line for line in script.lines))

    def test_failed_test_is_not_saved_unless_confirmed(self):
        code, text, _ = self._run(["1", "claude-sonnet-5-5", "y", "n"], "sk-ant-bad", lambda m: "401 invalid")
        self.assertEqual((code, text), (1, None))
        code, text, _ = self._run(["1", "claude-sonnet-5-5", "y", "y"], "sk-ant-bad", lambda m: "401 invalid")
        self.assertEqual(code, 0)
        self.assertIn(f"{LLM_ENV}=anthropic:claude-sonnet-5-5", text)

    def test_bad_choice_and_empty_key_save_nothing(self):
        self.assertEqual(self._run(["99"], "", lambda m: None)[:2], (2, None))
        self.assertEqual(self._run(["1", ""], "", lambda m: None)[:2], (2, None))
        self.assertEqual(self._run(["1", ""], "has space", lambda m: None)[:2], (2, None))

    def test_not_installed_provider_prints_install_command(self):
        groq_index = next(i for i, p in enumerate(PROVIDERS, start=1) if p.prefix == "groq")
        with mock.patch("gem_review.lit.setup_cli.is_installed", return_value=False):
            code, text, script = self._run([str(groq_index)], "", lambda m: None)
        self.assertEqual((code, text), (2, None))
        self.assertTrue(any("uv add 'pydantic-ai-slim[groq]'" in line for line in script.lines))


if __name__ == "__main__":
    unittest.main()
