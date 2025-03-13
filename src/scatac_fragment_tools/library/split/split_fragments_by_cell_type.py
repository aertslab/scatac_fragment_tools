from __future__ import annotations

import os
import shutil

from scatac_fragment_tools import _rust_scatac_fragment_tools



def _santize_string_for_filename(s: str) -> str:
    return s.replace(" ", "_").replace("/", "_")


def split_fragment_files_by_cell_type(
    sample_to_fragment_file: dict[str, str],
    path_to_temp_folder: str,
    path_to_output_folder: str,
    sample_to_cell_type_to_cell_barcodes: dict[str, dict[str, list]],
    chromsizes: dict[str, int],
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
    """
    # Check whether the same samples were provided in sample_to_fragment_file and
    # sample_to_cell_type_to_cell_barcodes.
    if set(sample_to_fragment_file.keys()) != set(
        sample_to_cell_type_to_cell_barcodes.keys()
    ):
        raise ValueError(
            "sample_to_fragment_file and sample_to_cell_type_to_cell_barcodes must have the same keys." + \
            f"\nsample_to_fragment_file has: {set(sample_to_fragment_file.keys())}" + \
            f"\nsample_to_cell_type_to_cell_barcodes has: {set(sample_to_cell_type_to_cell_barcodes.keys())}"
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

    cell_type_to_fragment_file_to_cell_barcode = {}
    for sample in sample_to_cell_type_to_cell_barcodes:
        fragment_file_path = sample_to_fragment_file[sample]
        for cell_type in sample_to_cell_type_to_cell_barcodes[sample]:
            if cell_type not in cell_type_to_fragment_file_to_cell_barcode:
                cell_type_to_fragment_file_to_cell_barcode[cell_type] = {
                    fragment_file_path: sample_to_cell_type_to_cell_barcodes[sample][cell_type]
                }
            else:
                cell_type_to_fragment_file_to_cell_barcode[cell_type] = {
                    **cell_type_to_fragment_file_to_cell_barcode[cell_type],
                    fragment_file_path: sample_to_cell_type_to_cell_barcodes[sample][cell_type]
                }

    _rust_scatac_fragment_tools.split_fragment_files_by_cell_type(
        fragment_file_paths=list(sample_to_fragment_file.values()),
        output_directory=path_to_output_folder,
        temp_directory=path_to_temp_folder + ("/" if not path_to_temp_folder.endswith("/") else ""),
        cell_type_to_fragment_file_to_cell_barcode=cell_type_to_fragment_file_to_cell_barcode,
        chromosomes=list(chromsizes.keys())
    )

    if clear_temp_folder:
        for cell_type in cell_type_to_fragment_file_to_cell_barcode:
            for chromosome in chromsizes.keys():
                os.remove(os.path.join(path_to_temp_folder, f"{cell_type}.{chromosome}.tsv.gz"))


    
    