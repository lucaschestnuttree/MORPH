#!/usr/bin/env python3
"""Prepare AlphaMissense tables already downloaded by the Wget wrapper."""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path
from typing import Sequence


REQUIRED_TABLES = (
    "AlphaMissense_hg38.tsv",
    "AlphaMissense_aa_substitutions.tsv",
)
DEFAULT_MISSENSE_DIR = Path("data/1-raw/missense_db")


def preprocess_file(file_path: Path) -> None:
    """Remove the three introductory lines from an uncompressed table."""
    temporary_path = file_path.with_name(file_path.name + ".tmp")
    try:
        with file_path.open("r") as source, temporary_path.open("w") as destination:
            for _ in range(3):
                if source.readline() == "":
                    raise ValueError(f"File has fewer than three introductory lines: {file_path}")
            shutil.copyfileobj(source, destination)
        temporary_path.replace(file_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def prepare_table(raw_dir: Path, table_name: str) -> bool:
    """Unzip and preprocess a downloaded table, or report that it is missing."""
    table_path = raw_dir / table_name
    compressed_path = raw_dir / f"{table_name}.gz"

    if table_path.is_file():
        print(f"[SKIP] {table_name} is already prepared.")
        return True
    if not compressed_path.is_file():
        print(f"[ERROR] Missing download: {compressed_path}")
        return False

    print(f"[UNZIP] {compressed_path.name}")
    try:
        with gzip.open(compressed_path, "rb") as source, table_path.open("wb") as destination:
            shutil.copyfileobj(source, destination)
        preprocess_file(table_path)
        compressed_path.unlink()
        print(f"[CLEANUP] Removed {compressed_path.name}")
        print(f"[READY] {table_path}")
        return True
    except (OSError, EOFError, ValueError) as error:
        table_path.unlink(missing_ok=True)
        print(f"[ERROR] Preparing {compressed_path}: {error}", file=sys.stderr)
        return False


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "base_path", nargs="?", type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Project root (defaults to the repository containing this script).",
    )
    parser.add_argument(
        "missense_dir", nargs="?", type=Path, default=DEFAULT_MISSENSE_DIR,
        help="Raw-data folder, relative to the project root unless absolute.",
    )
    args = parser.parse_args(argv)
    base_path = args.base_path.resolve()
    raw_dir = args.missense_dir if args.missense_dir.is_absolute() else base_path / args.missense_dir
    raw_dir = raw_dir.resolve()

    if not base_path.is_dir():
        print(f"[ERROR] Project root does not exist: {base_path}", file=sys.stderr)
        return 1
    raw_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Raw-data folder: {raw_dir}")
    failures = [name for name in REQUIRED_TABLES if not prepare_table(raw_dir, name)]
    if failures:
        print("[ERROR] Could not prepare: " + ", ".join(failures), file=sys.stderr)
        return 1

    print("[DONE] Required AlphaMissense tables are ready.")
    return 0


if __name__ == "__main__":
    sys.exit(main())