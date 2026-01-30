# Changelog

## v0.1.5

### Fixed
- **[6e3dbb9](https://github.com/aertslab/scatac_fragment_tools/commit/6e3dbb93ff9a459e8c920fae45f865d0c875148e)** Support reading CellRanger-ATAC 2.2.x and CellRanger-ARC 2.1.x fragment files.
- **[529a1c0](https://github.com/aertslab/scatac_fragment_tools/commit/529a1c0d0ccf1a214f1a27127e450416ce05c2a6)** Fix "--chrom_prefix" option of `scatac_fragment_tools bigwig`.

## v0.1.4

### Updated
- **[9e2ffaf](https://github.com/aertslab/scatac_fragment_tools/commit/9e2ffaf00d8cfe9148b7ae42603a96ab083d9f03)**: Add option to prefix cell barcode with sample id.
- **[5ddbda9](https://github.com/aertslab/scatac_fragment_tools/commit/5ddbda91a48cc7c4daec3b4d53b5dcf085bd179c)**: More memory efficient sorting of written out fragment files.

### Fixed

- **[3567a26](https://github.com/aertslab/scatac_fragment_tools/commit/3567a2632ba67a50730f5a83408cebbcc80c3a37)**: Don't panic when > 5 cols in fragment file.
- **[11c9e8f](https://github.com/aertslab/scatac_fragment_tools/commit/11c9e8ff6741a875a895f44cf116cde8099ee378)**: Update test code for newer polars version.


## v0.1.3

### Fixed

- **[f3603f02](https://github.com/aertslab/scatac_fragment_tools/commit/f3603f021a904a64b640e0c983366cb0e96a9c60)**: Throw error if `scatac_from_fragments split` can't write to fragments files, by updating rust_htslib as both `rust_htslib::bgzf:Reader` and`rust_htslib::bgzf:Writer` did not check if it could open a file successfully. Before this would silently fail later when trying to read those "written" files.


## v0.1.2

### Updated

- **[80db60fb](https://github.com/aertslab/scatac_fragment_tools/commit/80db60fbadeec22e08334f292134fa70e47c9019)**: Update Polars syntax to 1.0.0.
- **[27b3a059](https://github.com/aertslab/scatac_fragment_tools/commit/27b3a059c23462ebdb7240d10ce4be171cd67ce0)**: Update Polars `groupby` syntax to `group_by`.
- **[98a55b0d](https://github.com/aertslab/scatac_fragment_tools/commit/98a55b0de121883fca3660a9128871b67e3deead)**: Update pyo3 and itertools dependencies.

### Fixed

- **[80db60fb](https://github.com/aertslab/scatac_fragment_tools/commit/80db60fbadeec22e08334f292134fa70e47c9019)**: Fix running scatac_scatac_fragment_tools bigwig on Polars >= 1.0.0.
- **[d4d80c05](https://github.com/aertslab/scatac_fragment_tools/commit/d4d80c05840eceb362d87a036922300f61aaf9d1)**: Propagate verbose flag to fragments_to_coverage when specified in `scatac_scatac_fragment_tools bigwig`.
- **[65381a76](https://github.com/aertslab/scatac_fragment_tools/commit/65381a76adb07c3a657f325bb7bfcd06f636b244)**: Fix building for Linux ppc64le.


## v0.1.1

### Added

- **[004a265](https://github.com/aertslab/scatac_fragment_tools/commit/004a2654ecd5ed0a33be78f6fa5789c0a41deafb)**: Allow barcodes to map to multiple cell types while splitting fragments by cell type.
