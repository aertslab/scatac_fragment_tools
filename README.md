# scATAC Fragment Tools

[![PyPI](https://img.shields.io/pypi/v/scatac-fragment-tools?color=green)](https://pypi.org/project/scatac-fragment-tools/)
[![Docs](https://img.shields.io/badge/Docs-blue)](https://aertslab.github.io/scatac_fragment_tools/)

Tools for manipulating scATAC-seq fragment files.


## Installation

### From PyPI

```bash
# Install scatac-fragment-tools.
pip install scatac-fragment-tools

# Install scatac-fragment-tools with "pyBigWig" support (replacement for "pybigtools").
# ("pyBigWig" is slower than "pybigtools", so "pyBigWig" is only recommended if you
# would experience issues with bigWig file creation with "pybigtools".)
pip install scatac-fragment-tools[pybigwig]
```

### From source

```bash
# Install maturin.
pip install maturin

# Or install maturin with zig support (to build wheel that can run on older Linux distributions).
pip install maturin[zig]

# Clone scatac_fragments_tools git repo.
git clone https://github.com/aertslab/scatac_fragment_tools.git

cd scatac_fragment_tools

# Install in developer mode (aka `pip install -e .`)
maturin develop --release

# Or build a wheel and install afterwards with pip.
maturin build --release
# pip install rust/target/wheels/scatac_fragment_tools-${scatac_tools_version}-cp${python_version}-abi3-manylinux_${glibc_version}_x86_64.whl

# Or build a wheel that is more distributable (would still run on CentOS7).
maturin build --release --compatibility manylinux2014 --zig
# pip install rust/target/wheels/scatac_fragment_tools-${scatac_tools_version}-cp${python_version}-abi3-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
```


## Usage

Please visit the [documentation](https://aertslab.github.io/scatac_fragment_tools/)
