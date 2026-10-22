"""Which LLM provider runs the extraction and chain agents, and where its key comes from.

The model is a pydantic-ai model string, "<provider>:<model name>". It is chosen by
`--llm`, else `GEM_REVIEW_LLM` from platform/.env, else the Anthropic default. The key for the
chosen provider is read from the environment variable listed below (platform/.env is loaded
into the environment at CLI start without overriding variables already set).
"""

import importlib
import os
from dataclasses import dataclass
from pathlib import Path

from pydantic_ai.models.anthropic import AnthropicModelSettings
from pydantic_ai.settings import ModelSettings

from scripts.gem_annotate.model_layout import PLATFORM_ROOT

ENV_FILE = PLATFORM_ROOT / ".env"
LLM_ENV = "GEM_REVIEW_LLM"
DEFAULT_LLM = "anthropic:claude-opus-5-5"
MAX_TOKENS = 16000


@dataclass(frozen=True, slots=True)
class Provider:
    prefix: str              # pydantic-ai model-string prefix
    label: str
    key_vars: tuple[str, ...]  # any one of these supplies the key; empty = no key needed
    extra: str | None        # pydantic-ai-slim extra that installs the SDK
    module: str              # import name that proves the SDK is installed
    default_model: str
    key_url: str             # where the user gets a key


PROVIDERS: tuple[Provider, ...] = (
    Provider("anthropic", "Anthropic (Claude)", ("ANTHROPIC_API_KEY",), "anthropic", "anthropic",
             "claude-opus-5-5", "https://console.anthropic.com/settings/keys"),
    Provider("openai", "OpenAI (GPT)", ("OPENAI_API_KEY",), "openai", "openai",
             "gpt-5", "https://platform.openai.com/api-keys"),
    Provider("google-gla", "Google Gemini", ("GOOGLE_API_KEY", "GEMINI_API_KEY"), "google", "google.genai",
             "gemini-2.5-pro", "https://aistudio.google.com/apikey"),
    Provider("groq", "Groq", ("GROQ_API_KEY",), "groq", "groq",
             "llama-3.3-70b-versatile", "https://console.groq.com/keys"),
    Provider("mistral", "Mistral", ("MISTRAL_API_KEY",), "mistral", "mistralai",
             "mistral-large-latest", "https://console.mistral.ai/api-keys"),
    Provider("deepseek", "DeepSeek", ("DEEPSEEK_API_KEY",), "openai", "openai",
             "deepseek-chat", "https://platform.deepseek.com/api_keys"),
    Provider("openrouter", "OpenRouter (many models)", ("OPENROUTER_API_KEY",), "openai", "openai",
             "anthropic/claude-sonnet-4.5", "https://openrouter.ai/keys"),
    Provider("xai", "xAI (Grok)", ("XAI_API_KEY",), "xai", "xai_sdk",
             "grok-4", "https://console.x.ai"),
    Provider("ollama", "Ollama (local, no key)", (), "openai", "openai",
             "llama3.1", ""),
)


def provider_for(model_string: str) -> Provider | None:
    prefix = model_string.split(":", 1)[0]
    return next((p for p in PROVIDERS if p.prefix == prefix), None)


def is_installed(provider: Provider) -> bool:
    try:
        importlib.import_module(provider.module)
    except ImportError:
        return False
    return True


def key_present(provider: Provider) -> bool:
    return not provider.key_vars or any(os.environ.get(v) for v in provider.key_vars)


def load_env_file(path: Path = ENV_FILE) -> list[str]:
    """Put KEY=VALUE lines of path into os.environ (existing variables win). Returns the keys set."""
    if not path.exists():
        return []
    loaded = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip().removeprefix("export ").strip()
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
            value = value[1:-1]
        if key and key not in os.environ:
            os.environ[key] = value
            loaded.append(key)
    return loaded


def resolve_llm(flag: str | None) -> str:
    """--llm, else GEM_REVIEW_LLM, else the Anthropic default."""
    return flag or os.environ.get(LLM_ENV) or DEFAULT_LLM


def readiness(model_string: str) -> str | None:
    """None when the provider can run; otherwise one sentence saying what is missing."""
    if ":" not in model_string:
        return f"'{model_string}' is not a provider:model string; run `python -m gem_review.lit setup`."
    provider = provider_for(model_string)
    if provider is None:
        return None  # not in our table: let pydantic-ai decide
    if not is_installed(provider):
        return (f"{provider.label} is not installed; run `uv add 'pydantic-ai-slim[{provider.extra}]'` "
                "from platform/.")
    if not key_present(provider):
        return (f"No key for {provider.label}: set {' or '.join(provider.key_vars)} in platform/.env "
                "(run `python -m gem_review.lit setup`).")
    return None


def settings_for(model_string: str) -> ModelSettings:
    if model_string.startswith("anthropic:"):
        return AnthropicModelSettings(max_tokens=MAX_TOKENS, anthropic_effort="medium")
    return ModelSettings(max_tokens=MAX_TOKENS)


def cache_settings(model_string: str) -> ModelSettings | None:
    """Per-run settings that put a 5-minute cache breakpoint on the last static system block.

    The run passes the chunk block as a plain-string instruction, which pydantic-ai treats as
    static, so the breakpoint lands after the chunks (rules + chunks are the cached prefix; the
    question in the user message stays outside it). Other providers cache long prefixes on their
    own, so they get nothing extra here.
    """
    if model_string.startswith("anthropic:"):
        return AnthropicModelSettings(anthropic_cache_instructions="5m")
    return None


def write_env_values(path: Path, values: dict[str, str]) -> Path:
    """Set or replace KEY=VALUE lines, keep everything else, mode 600 (owner only)."""
    lines = path.read_text().splitlines() if path.exists() else []
    remaining = dict(values)
    out = []
    for line in lines:
        key = line.split("=", 1)[0].strip().removeprefix("export ").strip() if "=" in line else None
        if key in remaining:
            out.append(f"{key}={remaining.pop(key)}")
        else:
            out.append(line)
    out.extend(f"{key}={value}" for key, value in remaining.items())
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as handle:
        handle.write("\n".join(out) + "\n")
    os.chmod(tmp, 0o600)
    tmp.replace(path)
    return path
