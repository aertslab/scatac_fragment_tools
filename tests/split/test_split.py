import os
import pathlib

import polars as pl
from polars.testing import assert_frame_equal

TEST_DIRECTORY = pathlib.Path(__file__).parent.absolute()

def test_entrypoint():
    exit_status = os.system("scatac_fragment_tools split")
    assert exit_status == 0

def run_split_command(tmp_path, output_folder):
    os.system(f"cp {TEST_DIRECTORY}/a.fragments.tsv.gz {tmp_path}")
    os.system(f"cp {TEST_DIRECTORY}/a.fragments.tsv.gz.tbi {tmp_path}")
    os.system(f"cp {TEST_DIRECTORY}/b.fragments.tsv.gz {tmp_path}")
    os.system(f"cp {TEST_DIRECTORY}/b.fragments.tsv.gz.tbi {tmp_path}")
    os.system(f"cp {TEST_DIRECTORY}/sample_to_fragment.tsv {tmp_path}")
    os.system(f"cp {TEST_DIRECTORY}/cell_type_annotation.tsv {tmp_path}")
    os.system(f"cp {TEST_DIRECTORY}/hg38.chrom.sizes {tmp_path}")
    COMMAND = f"""cd {tmp_path} && \
    scatac_fragment_tools split \
        -f {tmp_path}/sample_to_fragment.tsv \
        -b {tmp_path}/cell_type_annotation.tsv \
        -c {tmp_path}/hg38.chrom.sizes \
        -o {output_folder} \
        -t {tmp_path} \
    """
    return os.system(COMMAND)

def test_split_command(tmp_path):
    output_folder = os.path.join(tmp_path, "output")
    os.makedirs(output_folder, exist_ok=True)
    exit_status = run_split_command(tmp_path, output_folder)
    assert exit_status == 0

    a_fragments = pl.read_csv(
        TEST_DIRECTORY.joinpath("a.fragments.tsv.gz"),
        separator = "\t",
        has_header = False
    )
    b_fragments = pl.read_csv(
        TEST_DIRECTORY.joinpath("b.fragments.tsv.gz"),
        separator = "\t",
        has_header = False
    )
    cell_annotations = pl.read_csv(
        TEST_DIRECTORY.joinpath("cell_type_annotation.tsv"),
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

