import os
import pathlib

import polars as pl
from polars.testing import assert_frame_equal

TEST_DIRECTORY = pathlib.Path(__file__).parent.absolute()

FILES_ALL_BARCODES_MAPPING_TO_SINGLE_TYPE = {
    "a.fragments": "a.fragments.tsv.gz",
    "a.fragment_index": "a.fragments.tsv.gz.tbi",
    "b.fragments": "b.fragments.tsv.gz",
    "b.fragment_index": "b.fragments.tsv.gz.tbi",
    "sample_to_fragment": "sample_to_fragment.tsv",
    "cell_type_annotation": "cell_type_annotation.tsv",
    "chrom_sizes": "hg38.chrom.sizes"
}

FILES_SOME_BARCODES_MAPPING_TO_MULTIPLE_TYPES = {
    "a.fragments": "a.fragments.tsv.gz",
    "a.fragment_index": "a.fragments.tsv.gz.tbi",
    "b.fragments": "b.fragments.tsv.gz",
    "b.fragment_index": "b.fragments.tsv.gz.tbi",
    "sample_to_fragment": "sample_to_fragment.tsv",
    "cell_type_annotation": "cell_type_annotation_one_bc_multiple_types.tsv",
    "chrom_sizes": "hg38.chrom.sizes"
}

def test_entrypoint():
    exit_status = os.system("scatac_fragment_tools split")
    assert exit_status == 0

def run_split_command(tmp_path, output_folder, file_dict):
    path_to_a_fragments = os.path.join(TEST_DIRECTORY, file_dict["a.fragments"])
    path_to_a_fragment_index = os.path.join(TEST_DIRECTORY, file_dict["a.fragment_index"])
    path_to_b_fragments = os.path.join(TEST_DIRECTORY, file_dict["b.fragments"])
    path_to_b_fragment_index = os.path.join(TEST_DIRECTORY, file_dict["b.fragment_index"])
    path_to_sample_to_fragment = os.path.join(TEST_DIRECTORY, file_dict["sample_to_fragment"])
    path_to_cell_type_annotation = os.path.join(TEST_DIRECTORY, file_dict["cell_type_annotation"])
    path_to_chrom_sizes = os.path.join(TEST_DIRECTORY, file_dict["chrom_sizes"])
    os.system(f"cp {path_to_a_fragments} {tmp_path}")
    os.system(f"cp {path_to_a_fragment_index} {tmp_path}")
    os.system(f"cp {path_to_b_fragments} {tmp_path}")
    os.system(f"cp {path_to_b_fragment_index} {tmp_path}")
    os.system(f"cp {path_to_sample_to_fragment} {tmp_path}")
    os.system(f"cp {path_to_cell_type_annotation} {tmp_path}")
    os.system(f"cp {path_to_chrom_sizes} {tmp_path}")

    COMMAND = f"""cd {tmp_path} && \
    scatac_fragment_tools split \
        -f {path_to_sample_to_fragment} \
        -b {path_to_cell_type_annotation} \
        -c {path_to_chrom_sizes} \
        -o {output_folder} \
        -t {tmp_path} \
    """
    return os.system(COMMAND)

def split_command_test_helper(tmp_path, file_dict):
    output_folder = os.path.join(tmp_path, "output")
    os.makedirs(output_folder, exist_ok=True)
    exit_status = run_split_command(tmp_path, output_folder, file_dict)
    assert exit_status == 0

    a_fragments = pl.read_csv(
        TEST_DIRECTORY.joinpath(file_dict["a.fragments"]),
        separator = "\t",
        has_header = False
    )
    b_fragments = pl.read_csv(
        TEST_DIRECTORY.joinpath(file_dict["b.fragments"]),
        separator = "\t",
        has_header = False
    )
    cell_annotations = pl.read_csv(
        TEST_DIRECTORY.joinpath(file_dict["cell_type_annotation"]),
        separator = "\t"
    )

    for row in cell_annotations \
            .select(pl.col("cell_type")) \
            .unique() \
            .iter_rows():
        cell_type = row[0]
        cell_type_bcs = cell_annotations \
            .filter(pl.col("cell_type") == cell_type) \
            .select(pl.col("cell_barcode"))
        fragments_cell_type = pl.concat(
            [
                a_fragments.filter(pl.col("column_4").is_in(cell_type_bcs)),
                b_fragments.filter(pl.col("column_4").is_in(cell_type_bcs))
            ]
        ).sort(by=["column_1", "column_2", "column_3", "column_4"])
        generated_fragments_cell_type = pl.read_csv(
            os.path.join(output_folder, f"{cell_type}.fragments.tsv.gz"),
            separator ="\t",
            has_header=False
        ).sort(by=["column_1", "column_2", "column_3", "column_4"])
        assert_frame_equal(
            fragments_cell_type,
            generated_fragments_cell_type
        )

def run_split_command_w_sample_id(tmp_path, output_folder, file_dict):
    path_to_a_fragments = os.path.join(TEST_DIRECTORY, file_dict["a.fragments"])
    path_to_a_fragment_index = os.path.join(TEST_DIRECTORY, file_dict["a.fragment_index"])
    path_to_b_fragments = os.path.join(TEST_DIRECTORY, file_dict["b.fragments"])
    path_to_b_fragment_index = os.path.join(TEST_DIRECTORY, file_dict["b.fragment_index"])
    path_to_sample_to_fragment = os.path.join(TEST_DIRECTORY, file_dict["sample_to_fragment"])
    path_to_cell_type_annotation = os.path.join(TEST_DIRECTORY, file_dict["cell_type_annotation"])
    path_to_chrom_sizes = os.path.join(TEST_DIRECTORY, file_dict["chrom_sizes"])
    os.system(f"cp {path_to_a_fragments} {tmp_path}")
    os.system(f"cp {path_to_a_fragment_index} {tmp_path}")
    os.system(f"cp {path_to_b_fragments} {tmp_path}")
    os.system(f"cp {path_to_b_fragment_index} {tmp_path}")
    os.system(f"cp {path_to_sample_to_fragment} {tmp_path}")
    os.system(f"cp {path_to_cell_type_annotation} {tmp_path}")
    os.system(f"cp {path_to_chrom_sizes} {tmp_path}")

    COMMAND = f"""cd {tmp_path} && \
    scatac_fragment_tools split \
        -f {path_to_sample_to_fragment} \
        -b {path_to_cell_type_annotation} \
        -c {path_to_chrom_sizes} \
        -o {output_folder} \
        -t {tmp_path} \
        --add_sample_id
    """
    return os.system(COMMAND)

def split_command_test_helper_w_sample_id(tmp_path, file_dict):
    output_folder = os.path.join(tmp_path, "output")
    os.makedirs(output_folder, exist_ok=True)
    exit_status = run_split_command_w_sample_id(tmp_path, output_folder, file_dict)
    assert exit_status == 0

    a_fragments = pl.read_csv(
        TEST_DIRECTORY.joinpath(file_dict["a.fragments"]),
        separator = "\t",
        has_header = False
    )
    b_fragments = pl.read_csv(
        TEST_DIRECTORY.joinpath(file_dict["b.fragments"]),
        separator = "\t",
        has_header = False
    )
    cell_annotations = pl.read_csv(
        TEST_DIRECTORY.joinpath(file_dict["cell_type_annotation"]),
        separator = "\t"
    )

    for (cell_type, ) in cell_annotations \
            .select([pl.col("cell_type")]) \
            .unique() \
            .iter_rows():
        cell_type_bcs = cell_annotations \
            .filter(pl.col("cell_type") == cell_type) \
            .select(pl.col("cell_barcode"))
        fragments_cell_type = pl.concat(
            [
                a_fragments.filter(
                    pl.col("column_4").is_in(cell_type_bcs)
                ).with_columns(
                    (pl.lit("A_") + pl.col("column_4")).alias("column_4")
                ),
                b_fragments.filter(pl.col("column_4").is_in(cell_type_bcs)
                ).with_columns(
                    (pl.lit("B_") + pl.col("column_4")).alias("column_4")
                )
            ]
        ).sort(by=["column_1", "column_2", "column_3", "column_4"])
        generated_fragments_cell_type = pl.read_csv(
            os.path.join(output_folder, f"{cell_type}.fragments.tsv.gz"),
            separator ="\t",
            has_header=False
        ).sort(by=["column_1", "column_2", "column_3", "column_4"])
        assert_frame_equal(
            fragments_cell_type,
            generated_fragments_cell_type
        )

def test_split_command_bc_single_type(tmp_path):
    split_command_test_helper(tmp_path, FILES_ALL_BARCODES_MAPPING_TO_SINGLE_TYPE)

def test_split_command_barcode_mapping_multiple_types(tmp_path):
    split_command_test_helper(tmp_path, FILES_SOME_BARCODES_MAPPING_TO_MULTIPLE_TYPES)

def test_split_command_w_sample_id_bc_single_type(tmp_path):
    split_command_test_helper_w_sample_id(tmp_path, FILES_ALL_BARCODES_MAPPING_TO_SINGLE_TYPE)

def test_split_command_w_sample_id_barcode_mapping_multiple_types(tmp_path):
    split_command_test_helper_w_sample_id(tmp_path, FILES_SOME_BARCODES_MAPPING_TO_MULTIPLE_TYPES)

