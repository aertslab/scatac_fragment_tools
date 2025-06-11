from __future__ import annotations

import argparse
import sys
from typing import Any, Callable, Dict, Optional

from rich_argparse import RichHelpFormatter

import scatac_fragment_tools
from scatac_fragment_tools.cli.commands import (
    command_fragment_to_bigwigs,
    command_split_fragments_by_cell_type,
)

# Constants
_VERSION = scatac_fragment_tools.__version__
_REQUIRED_ARGUMENTS_NAME = "Required arguments"
_OPTIONAL_ARGUMENTS_NAME = "Optional arguments"
_FORMATTER_CLASS = RichHelpFormatter


class SubparserBuilder:
    """
    Helper class for building subparsers.

    Automatically adds a required and an optional argument group
    and prints default values for optional arguments.
    Also, sets the default formatter class.
    """

    def __init__(
        self,
        name: str,
        subparsers: argparse._SubParsersAction,
        func: Optional[Callable] = None,
        **kwargs,
    ):
        self.name = name
        self.parser: argparse.ArgumentParser = subparsers.add_parser(
            name=name,
            formatter_class=_FORMATTER_CLASS,
            **kwargs,
        )
        self.parser.set_defaults(func=func)
        self.required_arguments = self.parser.add_argument_group(
            _REQUIRED_ARGUMENTS_NAME
        )
        self.optional_arguments = self.parser.add_argument_group(
            _OPTIONAL_ARGUMENTS_NAME
        )

    def add_required_argument(self, *args, **kwargs):
        self.required_arguments.add_argument(
            *args,
            required=True,
            **kwargs,
        )

    def add_optional_argument(self, *args, default: Any, help: str = "", **kwargs):
        self.optional_arguments.add_argument(
            *args,
            required=False,
            help=f"{help} Default: {default}",
            default=default,
            **kwargs,
        )

    def get_parser(self) -> Dict[str, argparse.ArgumentParser]:
        return {self.name: self.parser}


def add_fragments_to_bigwig_subparser(
    subparsers: argparse._SubParsersAction,
) -> Dict[str, argparse.ArgumentParser]:
    parser = SubparserBuilder(
        name="bigwig",
        subparsers=subparsers,
        func=command_fragment_to_bigwigs,
        description="Calculate genome coverage for fragments and write result to bigWig file.",
    )
    parser.add_required_argument(
        "-c",
        "--chrom",
        dest="chrom_sizes_filename",
        action="store",
        type=str,
        help="Filename with chromosome sizes (*.chrom.sizes, *.fa.fai).",
    )
    parser.add_required_argument(
        "-i",
        "--frag",
        dest="fragments_filename",
        action="store",
        type=str,
        help="Fragments TSV file for which to calculate genome coverage.",
    )
    parser.add_required_argument(
        "-o",
        "--bw",
        dest="bigwig_filename",
        action="store",
        type=str,
        help="BigWig filename to which the genome coverage data will be written.",
    )
    parser.add_optional_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store_true",
        default=False,
        help="Normalize genome coverage data by dividing by sequencing depth "
        "(number of fragments) multiplied by 1 million.",
    )
    parser.add_optional_argument(
        "-s",
        "--scaling",
        dest="scaling_factor",
        action="store",
        type=float,
        default=1.0,
        help="Scaling factor for genome coverage data. "
        "If normalization is enabled, scaling is applied afterwards.",
    )
    parser.add_optional_argument(
        "-x",
        "--cut-sites",
        dest="cut_sites",
        action="store_true",
        default=False,
        help="Use 1 bp Tn5 cut sites (start and end of each fragment) instead of whole "
        "fragment length for coverage calculation.",
    )
    parser.add_optional_argument(
        "-w",
        "--writer",
        dest="bigwig_writer",
        action="store",
        type=str.lower,
        choices=["pybigtools", "pybigwig"],
        default="pybigwig",
        help="Which bigWig writer implementation to use.",
    )
    parser.add_optional_argument(
        "--chrom-prefix",
        dest="chrom_prefix",
        action="store_true",
        default=False,
        help="Add chromosome prefix to each chromosome name found in the fragments file.",
    )
    parser.add_optional_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Whether to print progress.",
    )
    return parser.get_parser()


def add_split_fragments_by_cell_type_subparser(
    subparsers: argparse._SubParsersAction,
) -> Dict[str, argparse.ArgumentParser]:
    parser = SubparserBuilder(
        name="split",
        subparsers=subparsers,
        func=command_split_fragments_by_cell_type,
        description="Split fragment files by cell type.",
    )
    parser.add_required_argument(
        "-f",
        "--sample_fragments",
        dest="path_to_sample_to_fragment_definition",
        action="store",
        type=str,
        help="Path to a text file mapping sample names to fragment files.",
    )
    parser.add_required_argument(
        "-b",
        "--cell_type_barcodes",
        dest="path_to_cell_type_to_cell_barcode_definition",
        action="store",
        type=str,
        help="Path to a text file mapping samples to cell types and cell types to cell barcodes.",
    )
    parser.add_required_argument(
        "-c",
        "--chrom",
        dest="chrom_sizes_filename",
        action="store",
        type=str,
        help="Filename with chromosome sizes (*.chrom.sizes, *.fa.fai).",
    )
    parser.add_required_argument(
        "-o",
        "--output",
        dest="path_to_output_folder",
        action="store",
        type=str,
        help="Path to output folder.",
    )
    parser.add_optional_argument(
        "-t",
        "--temp",
        dest="path_to_temp_folder",
        action="store",
        type=str,
        default="/tmp",
        help="Path to temporary folder.",
    )
    parser.add_optional_argument(
        "-n",
        "--n_cpu",
        dest="n_cpu",
        action="store",
        type=int,
        default=1,
        help="Number of cores to use.",
    )
    parser.add_optional_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Whether to print progress.",
    )
    parser.add_optional_argument(
        "--clear_temp",
        dest="clear_temp_folder",
        action="store_true",
        default=False,
        help="Whether to clear the temporary folder.",
    )
    parser.add_optional_argument(
        "-s",
        "--sep",
        dest="separator",
        action="store",
        type=str,
        default="\t",
        help="Separator for text files.",
    )
    parser.add_optional_argument(
        "--sample_column",
        dest="sample_column_name",
        action="store",
        type=str,
        default="sample",
        help="Column name for the sample name",
    )
    parser.add_optional_argument(
        "--fragment_column",
        dest="path_to_fragment_file_column_name",
        action="store",
        type=str,
        default="path_to_fragment_file",
        help="Column name for the path to the fragment file",
    )
    parser.add_optional_argument(
        "--cell_type_column",
        dest="cell_type_column_name",
        action="store",
        type=str,
        default="cell_type",
        help="Column name for the cell type",
    )
    parser.add_optional_argument(
        "--cell_barcode_column",
        dest="cell_barcode_column_name",
        action="store",
        type=str,
        default="cell_barcode",
        help="Column name for the cell barcode",
    )
    return parser.get_parser()


_PARSERS_CREATOR_FUNCS = [
    add_fragments_to_bigwig_subparser,
    add_split_fragments_by_cell_type_subparser,
]


def main() -> int:
    main_parser = argparse.ArgumentParser(
        description=f"scATAC-fragment-tools (v{_VERSION}): Tools for processing scATAC-seq fragments.",
        formatter_class=RichHelpFormatter,
    )
    subparsers = main_parser.add_subparsers(
        title="subcommands",
    )
    _commands = {}
    for parser_creator_func in _PARSERS_CREATOR_FUNCS:
        _commands.update(parser_creator_func(subparsers))
    # When program is called without subcommand, print main help message.
    if len(sys.argv) == 1:
        main_parser.print_help()
    # When program is called with a subcommand, but no other arguments, print
    # subcommand's help message. When the subcommand is not recognized, print
    # main help message, unless the subcommand is "--help", in which case print.
    elif len(sys.argv) == 2:
        called_command = sys.argv[1]
        if called_command in _commands:
            _commands[called_command].print_help()
        elif called_command == "--help":
            main_parser.print_help()
        else:
            print(
                f"\033[1;31mUnrecognized command:\033[0m \033[1m{called_command}\033[0m"
            )
            main_parser.print_help()
    # When program is called with a subcommand and other arguments, parse arguments.
    else:
        args = main_parser.parse_args()
        if hasattr(args, "func"):
            args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
