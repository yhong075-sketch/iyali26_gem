"""`python -m gem_review.lit setup`: choose the provider and model, enter the key, test, save."""

import getpass
import os
from collections.abc import Callable
from pathlib import Path

from .providers import (ENV_FILE, LLM_ENV, PROVIDERS, Provider, is_installed, key_present, settings_for,
                        write_env_values)


def test_model(model_string: str) -> str | None:
    """One minimal call through pydantic-ai; None on success, else the error in one line."""
    from pydantic_ai import Agent

    try:
        agent = Agent(model_string, model_settings=settings_for(model_string))
        result = agent.run_sync("Reply with the single word OK.")
    except Exception as exc:  # any provider/network/auth error: shown to the user, never raised
        return f"{type(exc).__name__}: {str(exc)[:200]}"
    usage = result.usage
    return None if usage.requests else "no request was made"


def run_setup(*, env_path: Path = ENV_FILE, read: Callable[[str], str] = input,
              read_secret: Callable[[str], str] = getpass.getpass,
              tester: Callable[[str], str | None] = test_model, out: Callable[[str], None] = print) -> int:
    out("Choose the LLM provider for gem_review (keys are stored in platform/.env, owner-only):")
    for index, provider in enumerate(PROVIDERS, start=1):
        state = "installed" if is_installed(provider) else f"not installed: uv add 'pydantic-ai-slim[{provider.extra}]'"
        key = "no key needed" if not provider.key_vars else "key set" if key_present(provider) else "no key"
        out(f"  {index}. {provider.label:28} [{state}; {key}]")
    choice = read(f"Provider number [1-{len(PROVIDERS)}]: ").strip()
    if not choice.isdigit() or not 1 <= int(choice) <= len(PROVIDERS):
        out("No valid choice; nothing saved.")
        return 2
    provider: Provider = PROVIDERS[int(choice) - 1]
    if not is_installed(provider):
        out(f"{provider.label} is not installed. Run from platform/:\n  uv add 'pydantic-ai-slim[{provider.extra}]'\n"
            "then run setup again.")
        return 2
    model = read(f"Model name [{provider.default_model}]: ").strip() or provider.default_model
    model_string = f"{provider.prefix}:{model}"

    values = {LLM_ENV: model_string}
    if provider.key_vars:
        var = provider.key_vars[0]
        hint = f" (get one at {provider.key_url})" if provider.key_url else ""
        prompt = f"{var}{hint}; input is hidden"
        if key_present(provider):
            prompt += "; press Enter to keep the current key"
        secret = read_secret(prompt + ": ").strip()
        if secret:
            if any(ch.isspace() for ch in secret):
                out("The key contains whitespace; nothing saved.")
                return 2
            values[var] = secret
            os.environ[var] = secret  # so the test below and this process see it
        elif not key_present(provider):
            out("No key entered; nothing saved.")
            return 2

    if read("Test the model with one tiny call? [Y/n]: ").strip().lower() not in {"n", "no"}:
        problem = tester(model_string)
        if problem is None:
            out(f"OK: {model_string} answered.")
        else:
            out(f"Test failed: {problem}")
            if read("Save anyway? [y/N]: ").strip().lower() not in {"y", "yes"}:
                out("Nothing saved.")
                return 1

    write_env_values(env_path, values)
    saved = ", ".join(f"{k}={v}" if k == LLM_ENV else f"{k}=<hidden>" for k, v in values.items())
    out(f"Saved to {env_path}: {saved}")
    out("Every `python -m gem_review.lit` command now uses this model; override once with --llm provider:model.")
    return 0
