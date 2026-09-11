"""Process-local guards for the single-process build and static patch entrypoints."""

from contextlib import ExitStack, contextmanager
from functools import wraps
from inspect import signature
import socket
from unittest.mock import patch
import urllib.request

from cobra.util.solver import solvers
import requests


class ExecutionBlocked(BaseException):
    """Do not let annotation retry/diagnostic Exception handlers hide a leak."""


@contextmanager
def execution_limits(*, no_solve=False, allow_network=True):
    attempts = {"optimization": 0, "network": 0}

    def deny(kind):
        def blocked(*args, **kwargs):
            attempts[kind] += 1
            raise ExecutionBlocked(f"{kind} is disabled for this execution")

        return blocked

    # ponytail: process-local patches assume a single build thread; use process
    # isolation if these entrypoints are ever exposed by a concurrent server.
    with ExitStack() as stack:
        if no_solve:
            seen = set()
            for interface in solvers.values():
                for name in ("optimize", "_optimize"):
                    key = (interface.Model, name)
                    if key not in seen and hasattr(*key):
                        stack.enter_context(patch.object(*key, deny("optimization")))
                        seen.add(key)
            import swiglpk

            for name in ("glp_simplex", "glp_intopt", "glp_exact"):
                stack.enter_context(patch.object(swiglpk, name, deny("optimization")))
        if not allow_network:
            for obj, name in (
                (socket, "create_connection"),
                (socket, "getaddrinfo"),
                (socket.socket, "connect"),
                (socket.socket, "connect_ex"),
                (requests.sessions.Session, "request"),
                (urllib.request, "urlopen"),
            ):
                stack.enter_context(patch.object(obj, name, deny("network")))
        yield attempts


def guarded_execution(function):
    parameters = signature(function)

    @wraps(function)
    def guarded(*args, **kwargs):
        bound = parameters.bind(*args, **kwargs)
        bound.apply_defaults()
        candidate = bound.arguments.get("r608_curation_path") is not None
        with execution_limits(
            no_solve=candidate or bound.arguments.get("no_solve", False),
            allow_network=not candidate and bound.arguments.get("allow_network", True),
        ):
            return function(*args, **kwargs)

    return guarded
