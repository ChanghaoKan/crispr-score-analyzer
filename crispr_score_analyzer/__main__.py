"""Allow invoking the CLI via ``python -m crispr_score_analyzer``."""

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
