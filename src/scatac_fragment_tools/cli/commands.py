from __future__ import annotations

from typing import Dict


def command_fragment_to_bigwigs(args):
    import os

    # Check arguments before doing anything else.
    if not os.path.exists(args.chrom_sizes_filename):
        raise FileNotFoundError(
            f"Chromosome sizes file not found: {args.chrom_sizes_filename}"
        )
    if not os.path.exists(args.fragments_filename):
        raise FileNotFoundError(f"Fragments file not found: {args.fragments_filename}")

    # Import inside function to avoid unnecessary imports when calling other commands.
    import polars as pl

    from scatac_fragment_tools.library.bigwig.fragments_to_bigwig import (
        fragments_to_bw,
        get_chromosome_sizes,
        read_fragments_to_polars_df,
    )

    chrom_sizes = get_chromosome_sizes(args.chrom_sizes_filename)
    fragments_df = read_fragments_to_polars_df(args.fragments_filename)

    if args.chrom_prefix:
        fragments_df = fragments_df.with_columns(
            (pl.lit(f"{args.chrom_prefix}_") + pl.col("Chromosome").cast(pl.Utf8))
            .cast(pl.Categorical)
            .alias("Chromosome")
        )

    fragments_to_bw(
        fragments_df,
        chrom_sizes,
        args.bigwig_filename,
        args.normalize,
        args.scaling_factor,
        args.cut_sites,
        args.bigwig_writer,
        args.verbose,
    )


def command_split_fragments_by_cell_type(args):
    """
    Split fragment files by cell type.

    Parameters
    ----------
    args: Namespace
        Command line arguments.

    Arguments
    ---------
    args.path_to_sample_to_fragment_definition: str
        Path to a text file
    args.path_to_cell_type_to_cell_barcode_definition: str
        Path to a text file
    args.path_to_temp_folder: str
        Path to temporary folder.
    args.path_to_output_folder: str
        Path to output folder.
    args.chrom_sizes_filename: str
        Path to text file with chromosome sizes.
    args.n_cpu: int
        Number of cores to use.
    args.verbose: bool
        Whether to print progress.
    args.clear_temp_folder: bool
        Whether to clear the temporary folder.
    args.separator: str
        Separator for text files.
    args.sample_column_name: str
        Column name for the sample name
    args.path_to_fragment_file_column_name: str
        Column name for the path to the fragment file
    args.cell_type_column_name: str
        Column name for the cell type
    args.cell_barcode_column_name: str
        Column name for the cell barcode
    """
    # Check arguments before doing anything else.
    import os

    if not os.path.exists(args.path_to_sample_to_fragment_definition):
        raise FileNotFoundError(
            f"Sample to fragment definition file not found: {args.path_to_sample_to_fragment_definition}"
        )
    if not os.path.exists(args.path_to_cell_type_to_cell_barcode_definition):
        raise FileNotFoundError(
            f"Cell type to cell barcode definition file not found: {args.path_to_cell_type_to_cell_barcode_definition}"
        )
    if not os.path.exists(args.chrom_sizes_filename):
        raise FileNotFoundError(
            f"Chromosome sizes file not found: {args.chrom_sizes_filename}"
        )

    # check headers of definition files, by reading the first line
    with open(args.path_to_sample_to_fragment_definition) as f:
        sample_to_fragment_header = f.readline().strip().split(args.separator)
        if not (
            args.sample_column_name in sample_to_fragment_header
            and args.path_to_fragment_file_column_name in sample_to_fragment_header
        ):
            raise KeyError(
                f"Sample to fragment definition file must have columns {args.sample_column_name} and {args.path_to_fragment_file_column_name}"
            )
    with open(args.path_to_cell_type_to_cell_barcode_definition) as f:
        cell_type_to_cell_barcode_header = f.readline().strip().split(args.separator)
        if not (
            args.cell_type_column_name in cell_type_to_cell_barcode_header
            and args.cell_barcode_column_name in cell_type_to_cell_barcode_header
            and args.sample_column_name in cell_type_to_cell_barcode_header
        ):
            raise KeyError(
                f"Cell type to cell barcode definition file must have columns {args.cell_type_column_name}, {args.cell_barcode_column_name} and {args.sample_column_name}"
            )

    import polars as pl

    from scatac_fragment_tools.library.split.split_fragments_by_cell_type import (
        split_fragment_files_by_cell_type,
    )

    # Read sample to fragment file definition
    # and create a dictionary mapping sample names to fragment files.
    sample_to_fragment_file: Dict[str, str] = {}
    d_sample_to_fragment_definition = pl.read_csv(
        args.path_to_sample_to_fragment_definition,
        separator=args.separator,
    ).to_dict()
    for sample, fragment_file_path in zip(
        d_sample_to_fragment_definition[args.sample_column_name],
        d_sample_to_fragment_definition[args.path_to_fragment_file_column_name],
    ):
        if sample in sample_to_fragment_file:
            raise ValueError(
                f"Duplicate sample name: {sample} in {args.path_to_sample_to_fragment_definition}"
            )
        sample_to_fragment_file[sample] = fragment_file_path

    # Read cell type to cell barcode definition
    # and create a dictionary mapping sample names to cell type to list of cell barcodes.
    sample_to_cell_type_to_cell_barcodes: Dict[str, Dict[str, list]] = {}
    d_cell_type_to_cell_barcode_definition = (
        pl.scan_csv(
            args.path_to_cell_type_to_cell_barcode_definition, separator=args.separator
        )
        .group_by([args.sample_column_name, args.cell_type_column_name])
        .agg(pl.col(args.cell_barcode_column_name))
        .collect()
        .to_dict()
    )
    for sample, cell_type, cell_barcodes in zip(
        d_cell_type_to_cell_barcode_definition[args.sample_column_name],
        d_cell_type_to_cell_barcode_definition[args.cell_type_column_name],
        d_cell_type_to_cell_barcode_definition[args.cell_barcode_column_name],
    ):
        if sample not in sample_to_cell_type_to_cell_barcodes:
            sample_to_cell_type_to_cell_barcodes[sample] = {}
        if cell_type in sample_to_cell_type_to_cell_barcodes[sample]:
            raise ValueError(
                f"Duplicates in {args.path_to_cell_type_to_cell_barcode_definition}, for sample {sample} and cell type {cell_type}"
            )
        sample_to_cell_type_to_cell_barcodes[sample][cell_type] = cell_barcodes

    # Read chromosome sizes
    # and create a dictionary mapping chromosome names to chromosome sizes.
    chromsizes: Dict[str, int] = {}
    pl_chromsizes = pl.read_csv(
        args.chrom_sizes_filename,
        separator="\t",
        has_header=False,
    )
    for chromosome, size in pl_chromsizes.iter_rows():
        if chromosome in chromsizes:
            raise ValueError(
                f"Duplicates in {args.chrom_sizes_filename}, for chromosome {chromosome}"
            )
        chromsizes[chromosome] = size

    split_fragment_files_by_cell_type(
        sample_to_fragment_file=sample_to_fragment_file,
        path_to_temp_folder=args.path_to_temp_folder,
        path_to_output_folder=args.path_to_output_folder,
        sample_to_cell_type_to_cell_barcodes=sample_to_cell_type_to_cell_barcodes,
        chromsizes=chromsizes,
        n_cpu=args.n_cpu,
        verbose=args.verbose,
        clear_temp_folder=args.clear_temp_folder,
    )
