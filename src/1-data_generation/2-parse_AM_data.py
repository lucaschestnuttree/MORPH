#!/usr/bin/env python3
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import polars as pl


@dataclass
class SplitContext:
    """
    Context object holding paths and naming info for the split operation.
    """
    base_path: Path
    tsv_path: Path
    suffix: str        # e.g. "_is_afile.tsv"
    stem_name: str     # e.g. "is_afile"
    out_dir: Path      # e.g. base/data/1-raw/missense_jul/is_afile


def extract_suffix_and_stem(tsv_path: Path) -> Tuple[str, str]:
    """
    Given a TSV path like 'this_is_afile.tsv', return:

        suffix:    '_is_afile.tsv'
        stem_name: 'is_afile'

    If there is no underscore in the filename, we fall back to:

        suffix:    '.tsv'
        stem_name: <name_without_extension>
    """
    name = tsv_path.name  # e.g. 'this_is_afile.tsv'

    if not name.endswith(".tsv"):
        raise ValueError(f"Input file does not end with .tsv: {name}")

    first_underscore = name.find("_")

    if first_underscore == -1:
        stem_name = tsv_path.stem
        suffix = ".tsv"
        return suffix, stem_name

    suffix = name[first_underscore:]          # "_is_afile.tsv"
    stem_name = suffix.lstrip("_")[:-4]       # "is_afile"
    return suffix, stem_name


def create_context(base_path_str: str, tsv_path_str: str) -> SplitContext:
    """
    Factory-style constructor for SplitContext.
    Normalises paths, derives suffix/stem, and creates output directory.
    """
    base_path = Path(base_path_str).resolve()
    tsv_path = Path(tsv_path_str).resolve()

    if not tsv_path.is_file():
        raise FileNotFoundError(f"Input TSV file not found: {tsv_path}")

    suffix, stem_name = extract_suffix_and_stem(tsv_path)

    out_dir = base_path / "data" / "2-processed" / "missense_jul" / stem_name
    out_dir.mkdir(parents=True, exist_ok=True)

    return SplitContext(
        base_path=base_path,
        tsv_path=tsv_path,
        suffix=suffix,
        stem_name=stem_name,
        out_dir=out_dir,
    )


def build_lazy_frame(context: SplitContext) -> pl.LazyFrame:
    """
    Build a Polars LazyFrame for the input TSV.

    This uses lazy scanning so Polars can optimise and parallelise parsing.
    """
    return pl.scan_csv(
        context.tsv_path,
        separator="\t",      # <-- use 'separator' instead of 'sep'
        has_header=True,
        infer_schema_length=1000,  # tweak if needed
    )


def split_by_uniprot(
    context: SplitContext,
    lazy_frame: pl.LazyFrame,
    uniprot_column: str = "uniprot_id",
) -> None:
    """
    Split the TSV into separate per-uniprot files.

    For each distinct value in `uniprot_column`, we create a TSV file named:

        <uniprot_value><suffix>

    in the directory:

        context.out_dir

    Each file contains all columns, but only rows with that uniprot_id.
    """
    df = lazy_frame.collect()

    if uniprot_column not in df.columns:
        cols = ", ".join(df.columns)
        msg = (
            f"Column '{uniprot_column}' not found in TSV. "
            f"Available columns: {cols}"
        )
        raise KeyError(msg)

    grouped = df.group_by(uniprot_column, maintain_order=True)

    for group_key, group_df in grouped:
        # Polars gives a tuple of key values, even for a single key
        if isinstance(group_key, tuple):
            uniprot_str = str(group_key[0])
        else:
            uniprot_str = str(group_key)

        file_name = f"{uniprot_str}{context.suffix}"
        output_path = context.out_dir / file_name

        group_df.write_csv(
            output_path,
            separator="\t",
            include_header=True,
        )


if __name__ == "__main__":
    # -------------------------------------------
    # Argument handling HERE ONLY (nowhere else)
    # -------------------------------------------

    if len(sys.argv) != 3:
        prog = Path(sys.argv[0]).name
        msg = (
            f"Usage:\n"
            f"  {prog} BASE_PATH TSV_PATH\n\n"
            f"  BASE_PATH  Root directory for output "
            f"(e.g. /cluster/project/root)\n"
            f"  TSV_PATH   Path to the large input TSV file\n"
        )
        print(msg, file=sys.stderr)
        sys.exit(1)

    base_path_str = sys.argv[1]
    tsv_path_str = sys.argv[2]

    # -------------------------------------------
    # Normal workflow
    # -------------------------------------------

    context = create_context(base_path_str, tsv_path_str)
    lazy_frame = build_lazy_frame(context)
    split_by_uniprot(context, lazy_frame)
