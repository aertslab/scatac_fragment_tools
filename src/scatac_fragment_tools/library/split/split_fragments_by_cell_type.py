from __future__ import annotations

import os
from typing import Dict, List

import joblib

from scatac_fragment_tools import _rust_scatac_fragment_tools

NUMBER_OF_WRITER_THREADS = 5


def _santize_string_for_filename(s: str) -> str:
    return s.replace(" ", "_").replace("/", "_")


def split_fragment_files_by_cell_type(
    sample_to_fragment_file: Dict[str, str],
    path_to_temp_folder: str,
    path_to_output_folder: str,
    sample_to_cell_type_to_cell_barcodes: Dict[str, Dict[str, list]],
    chromsizes: Dict[str, int],
    n_cpu: int = 1,
    verbose: bool = False,
    clear_temp_folder: bool = False,
):
    """
    Split fragment files by cell type.

    Parameters
    ----------
    sample_to_fragment_file : Dict[str, str]
        Dictionary mapping sample names to fragment files.
    path_to_temp_folder : str
        Path to temporary folder, used for writing fragment files
        per cell type split by sample.
    path_to_output_folder : str
        Path to output folder, used for writing fragment files
        per cell type (merged across samples).
    sample_to_cell_type_to_cell_barcodes : Dict[str, Dict[str, list]]
        Dictionary mapping sample names to cell type to list of cell barcodes.
    chromsizes : Dict[str, int]
        Dictionary mapping chromosome names to chromosome sizes.
    n_cpu : int, optional
        Number of corse to use. The default is 1.
    verbose : bool, optional
        Whether to print progress. The default is False.
    clear_temp_folder : bool, optional
        Whether to clear the temporary folder. The default is False.
    """
    # Check whether the same samples were provided in sample_to_fragment_file and
    # sample_to_cell_type_to_cell_barcodes.
    if set(sample_to_fragment_file.keys()) != set(
        sample_to_cell_type_to_cell_barcodes.keys()
    ):
        raise ValueError(
            "sample_to_fragment_file and sample_to_cell_type_to_cell_barcodes must have the same keys."
        )

    # Create tmp folder if it does not exist.
    if not os.path.exists(path_to_temp_folder):
        if verbose:
            print(f"Creating {path_to_temp_folder}")
        os.makedirs(path_to_temp_folder)

    # Create output folder if it does not exists.
    if not os.path.exists(path_to_output_folder):
        if verbose:
            print(f"Creating {path_to_output_folder}")
        os.makedirs(path_to_output_folder)

    # Create a folder for each sample.
    for sample in sample_to_fragment_file:
        os.makedirs(os.path.join(path_to_temp_folder, sample), exist_ok=True)

    # Split fragment files by cell barcode, in parallel.
    if verbose:
        print("Splitting fragments ...")
    joblib.Parallel(n_jobs=n_cpu)(
        joblib.delayed(_rust_scatac_fragment_tools.split_fragments_by_cell_barcode)(
            path_to_fragments=sample_to_fragment_file[sample],
            path_to_output_folder=os.path.join(path_to_temp_folder, sample),
            cell_type_to_cell_barcodes=sample_to_cell_type_to_cell_barcodes[sample],
            chromsizes=chromsizes,
            verbose=verbose,
        )
        for sample in sample_to_cell_type_to_cell_barcodes
    )

    # Check whether all files were create successfully and
    # create a dictionary mapping cell types to fragment files.
    cell_type_to_fragment_files: Dict[str, List[str]] = {}
    for sample in sample_to_cell_type_to_cell_barcodes:
        for cell_type in sample_to_cell_type_to_cell_barcodes[sample]:
            cell_type_sanitized = _santize_string_for_filename(cell_type)
            path_to_fragment_file = os.path.join(
                path_to_temp_folder, sample, f"{cell_type_sanitized}.fragments.tsv.gz"
            )
            if not os.path.exists(path_to_fragment_file):
                raise ValueError(
                    f"Fragment file {path_to_fragment_file} does not exist."
                )
            if cell_type_sanitized not in cell_type_to_fragment_files:
                cell_type_to_fragment_files[cell_type_sanitized] = []
            cell_type_to_fragment_files[cell_type_sanitized].append(
                path_to_fragment_file
            )

    # Merge fragment files by cell type, in parallel.
    if verbose:
        print("Merging fragments ...")
    joblib.Parallel(n_jobs=n_cpu)(
        joblib.delayed(_rust_scatac_fragment_tools.merge_fragment_files)(
            path_to_fragment_files=cell_type_to_fragment_files[cell_type],
            path_to_output_file=os.path.join(
                path_to_output_folder, f"{cell_type}.fragments.tsv.gz"
            ),
            number_of_threads=NUMBER_OF_WRITER_THREADS,
            verbose=verbose,
        )
        for cell_type in cell_type_to_fragment_files
    )

    # Check whether all files were create successfully.
    for cell_type in cell_type_to_fragment_files:
        cell_type_sanitized = _santize_string_for_filename(cell_type)
        path_to_fragment_file = os.path.join(
            path_to_output_folder, f"{cell_type_sanitized}.fragments.tsv.gz"
        )
        if not os.path.exists(path_to_fragment_file):
            Warning(f"Fragment file {path_to_fragment_file} does not exist.")

    # Clear temporary folder.
    if clear_temp_folder:
        if verbose:
            print("Clearing temporary folder ...")
        for cell_type in cell_type_to_fragment_files:
            for fragment_file in cell_type_to_fragment_files[cell_type]:
                if fragment_file.startswith(
                    path_to_temp_folder
                ) and fragment_file.endswith(".fragments.tsv.gz"):
                    if verbose:
                        print(f"Removing {fragment_file}")
                    os.remove(fragment_file)
