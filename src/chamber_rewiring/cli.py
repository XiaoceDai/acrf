"""Console entry points: run bundled scripts with correct cwd and PYTHONPATH."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from chamber_rewiring.paths import project_root


def _env_with_root(root: Path) -> dict[str, str]:
    env = os.environ.copy()
    root_s = str(root)
    prev = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = root_s if not prev else f"{root_s}{os.pathsep}{prev}"
    return env


def cmd_release_bundle(_: argparse.Namespace) -> None:
    root = project_root()
    script = root / "scripts" / "innovation_release_bundle.py"
    subprocess.run(
        [sys.executable, str(script)],
        cwd=root,
        env=_env_with_root(root),
        check=True,
    )


def cmd_run(ns: argparse.Namespace) -> None:
    root = project_root()
    script = Path(ns.script)
    if not script.is_absolute():
        script = root / script
    if not script.is_file():
        raise SystemExit(f"not a file: {script}")
    subprocess.run(
        [sys.executable, str(script), *ns.args],
        cwd=root,
        env=_env_with_root(root),
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="chamber-rewiring",
        description="Project P02 pipeline helpers (install with pip install -e .).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_bundle = sub.add_parser(
        "release-bundle",
        help="Run scripts/innovation_release_bundle.py (NCP v3 + Fig2 pack + manifest).",
    )
    p_bundle.set_defaults(func=cmd_release_bundle)

    p_run = sub.add_parser(
        "run",
        help="Run any project script under repo root (sets PYTHONPATH to include root).",
    )
    p_run.add_argument(
        "script",
        help="Path relative to repo root, e.g. scripts/fig6b_v3_anchor_matrix.py",
    )
    p_run.add_argument(
        "args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to the script.",
    )
    p_run.set_defaults(func=cmd_run)

    ns = parser.parse_args()
    ns.func(ns)


if __name__ == "__main__":
    main()
