import scatac_fragment_tools
from scatac_fragment_tools.cli.commands import command_fragment_to_bigwigs
import argparse
from rich_argparse import RichHelpFormatter
import sys
from typing import Dict

# Constants
_VERSION = scatac_fragment_tools.__version__


def add_fragments_to_bigwig_subparser(
    subparsers: argparse._SubParsersAction,
) -> Dict[str, argparse.ArgumentParser]:
    _NAME = "bigwig"
    parser = subparsers.add_parser(
        name=_NAME,
        add_help=True,
        description="Calculate genome coverage for fragments and write result to bigWig file.",
        formatter_class=RichHelpFormatter,
    )
    parser.set_defaults(func=command_fragment_to_bigwigs)
    required_arguments = parser.add_argument_group("Required arguments")
    optional_arguments = parser.add_argument_group("Optional arguments")
    required_arguments.add_argument(
        "-c",
        "--chrom",
        dest="chrom_sizes_filename",
        action="store",
        type=str,
        required=True,
        help="Filename with chromosome sizes (*.chrom.sizes, *.fa.fai).",
    )
    required_arguments.add_argument(
        "-i",
        "--frag",
        dest="fragments_filename",
        action="store",
        type=str,
        required=True,
        help="Fragments TSV file for which to calculate genome coverage.",
    )
    required_arguments.add_argument(
        "-o",
        "--bw",
        dest="bigwig_filename",
        action="store",
        type=str,
        required=True,
        help="BigWig filename to which the genome coverage data will be written.",
    )
    optional_arguments.add_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store",
        type=str,
        required=False,
        choices={"yes", "no"},
        default="yes",
        help='Normalize genome coverage data. Default: "yes".',
    )
    optional_arguments.add_argument(
        "-s",
        "--scaling",
        dest="scaling_factor",
        action="store",
        type=float,
        required=False,
        default=1.0,
        help='Scaling factor for genome coverage data. Default: "1.0".',
    )
    optional_arguments.add_argument(
        "--chrom-prefix",
        dest="chrom_prefix",
        action="store",
        type=str,
        required=False,
        help="Add chromosome prefix to each chromosome name found in the fragments file.",
    )
    return {_NAME: parser}


_PARSERS_CREATOR_FUNCS = [
    add_fragments_to_bigwig_subparser,
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
    # When program is called without subcommand, print main help message
    if len(sys.argv) == 1:
        main_parser.print_help()
    # When program is called with a subcommand, but no other arguments, print subcommand help message
    # If the subcommand is not recognized, print main help message
    elif len(sys.argv) == 2:
        called_command = sys.argv[1]
        if called_command in _commands:
            _commands[called_command].print_help()
        else:
            print(
                f"\033[1;31mUnrecognized command:\033[0m \033[1m{called_command}\033[0m"
            )
            main_parser.print_help()
    # When program is called with a subcommand and other arguments, parse arguments
    else:
        args = main_parser.parse_args()
        if hasattr(args, "func"):
            args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
