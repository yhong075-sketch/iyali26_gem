"""iYali26 reference build entry point."""

def main(*args, **kwargs):
    from .cli import main as run_cli
    return run_cli(*args, **kwargs)

__all__ = ["main"]
