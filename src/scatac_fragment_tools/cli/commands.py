def command_fragment_to_bigwigs(args):
    import os

    # Check arguments before doing anything else.
    if not os.path.exists(args.chrom_sizes_filename):
        raise FileNotFoundError(
            f"Chromosome sizes file not found: {args.chrom_sizes_filename}"
        )
    if not os.path.exists(args.fragments_filename):
        raise FileNotFoundError(f"Fragments file not found: {args.fragments_filename}")

    # import inside function to avoid unnecessary imports when calling other commands.
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
    )
