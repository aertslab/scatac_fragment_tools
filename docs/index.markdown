---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
nav_order: 1
---

# scATAC-seq Fragment Tools

Tools for manipulating scATAC-seq fragment files.

## Installation

{: .note}
> Requires python version **>=3.8**.

### From pypi:

{% highlight bash %}
# Install scatac-fragment-tools.
pip install scatac-fragment-tools

# Install scatac-fragment-tools with "pyBigWig" support (replacement for "pybigtools").
# ("pyBigWig" is slower than "pybigtools", so "pyBigWig" is only recommended if you
# would experience issues with bigWig file creation with "pybigtools".)
pip install scatac-fragment-tools[pybigwig]

{% endhighlight %}

### From source:

{% highlight bash %}
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
#pip install rust/target/wheels/scatac_fragment_tools-${scatac_tools_version}-cp${python_version}-abi3-manylinux_2_17_x86_64.manylinux2014_x86_64.whl

{% endhighlight %}

## Testing installation

{% highlight bash %}
> scatac_fragment_tools
Usage: scatac_fragment_tools [-h] {bigwig,split} ...

scATAC-fragment-tools (v0.1.5): Tools for processing scATAC-seq fragments.

Options:
  -h, --help      show this help message and exit

Subcommands:
  {bigwig,split}

{% endhighlight %}
