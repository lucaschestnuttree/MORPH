#!/usr/bin/env python3
"""Split prepared AlphaMissense TSV tables into one file per UniProt protein.

The raw input tables must already exist in data/1-raw/missense_db/. Run
the server_wrappers/_0-fetch.sh wrapper with sbatch first when either table is missing.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import polars as pl


DEFAULT_TABLES = (
    "AlphaMissense_aa_substitutions.tsv",
    "AlphaMissense_hg38.tsv",
)


@dataclass
class SplitContext:
    base_path: Path
    tsv_path: Path
    suffix: str
    stem_name: str
    out_dir: Path


def extract_suffix_and_stem(tsv_path: Path) -> tuple[str, str]:
    """Keep the original per-protein naming convention derived from the input name."""
    name = tsv_path.name
    if not name.endswith(".tsv"):
        raise ValueError(f"Input file does not end with .tsv: {name}")
    first_underscore = name.find("_")
    if first_underscore == -1:
        return ".tsv", tsv_path.stem
    suffix = name[first_underscore:]
    return suffix, suffix.lstrip("_")[:-4]


def create_context(base_path: Path, tsv_path: Path) -> SplitContext:
    if not tsv_path.is_file():
        raise FileNotFoundError(f"Input TSV file not found: {tsv_path}")
    suffix, stem_name = extract_suffix_and_stem(tsv_path)
    out_dir = base_path / "data" / "2-processed" / "missense_jul" / stem_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return SplitContext(base_path, tsv_path, suffix, stem_name, out_dir)


def build_lazy_frame(context: SplitContext) -> pl.LazyFrame:
    return pl.scan_csv(
        context.tsv_path,
        separator="\t",
        has_header=True,
        infer_schema_length=1000,
    )


def split_by_uniprot(
    context: SplitContext,
    lazy_frame: pl.LazyFrame,
    uniprot_column: str = "uniprot_id",
) -> int:
    """Write all columns and rows for each protein, preserving the original split."""
    dataframe = lazy_frame.collect()
    if uniprot_column not in dataframe.columns:
        columns = ", ".join(dataframe.columns)
        raise KeyError(
            f"Column '{uniprot_column}' not found in {context.tsv_path}. "
            f"Available columns: {columns}"
        )

    count = 0
    for group_key, group_dataframe in dataframe.group_by(uniprot_column, maintain_order=True):
        uniprot_id = str(group_key[0] if isinstance(group_key, tuple) else group_key)
        output_path = context.out_dir / f"{uniprot_id}{context.suffix}"
        group_dataframe.write_csv(output_path, separator="\t", include_header=True)
        count += 1
    return count



def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "base_path", nargs="?", default=Path(__file__).resolve().parents[2],
        type=Path, metavar="BASE_PATH",
        help="Project root; defaults to the repository containing this script.",
    )
    parser.add_argument(
        "--missense-dir", default="data/1-raw/missense_db", metavar="PATH",
        help="Raw-data directory, relative to BASE_PATH unless absolute.",
    )
    parser.add_argument(
        "--tsv", action="append", metavar="PATH",
        help="Split this TSV instead of the two default tables; repeat for several. "
             "Relative paths resolve against BASE_PATH.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_argument_parser().parse_args(argv)
    base_path = args.base_path.resolve()
    raw_dir = Path(args.missense_dir)
    if not raw_dir.is_absolute():
        raw_dir = base_path / raw_dir
    raw_dir = raw_dir.resolve()
    tsv_paths = (
        [Path(path) if Path(path).is_absolute() else base_path / path for path in args.tsv]
        if args.tsv else [raw_dir / name for name in DEFAULT_TABLES]
    )
    try:
        if not base_path.is_dir():
            raise FileNotFoundError(f"Project root not found: {base_path}")
        missing = [str(path) for path in tsv_paths if not path.is_file()]
        if missing:
            raise FileNotFoundError(
                "Input TSV file(s) not found: " + ", ".join(missing)
                + ". Submit server_wrappers/_0-fetch.sh with sbatch first if these are the raw AlphaMissense tables."
            )
        print(f"[INFO] Project root: {base_path}")
        print(f"[INFO] Raw-data folder: {raw_dir}")
        split_results = []
        for tsv_path in tsv_paths:
            context = create_context(base_path, tsv_path)
            print(f"[RUN] Splitting {tsv_path.name} by uniprot_id")
            count = split_by_uniprot(context, build_lazy_frame(context))
            split_results.append((context, count))
            print(f"[DONE] Wrote {count:,} per-protein TSVs to {context.out_dir}")
    except (OSError, ValueError, KeyError, EOFError, pl.exceptions.PolarsError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1

    print("\n[SUMMARY] Split TSV files written in this run:")
    for context, count in split_results:
        output_folder = context.out_dir.relative_to(base_path).as_posix()
        print(f"[SUMMARY] {context.stem_name}: {count:,} files -> {output_folder}/")
    total = sum(count for _, count in split_results)
    print(f"[SUMMARY] Total: {total:,} files")
    print("[DONE] AlphaMissense table splitting completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
