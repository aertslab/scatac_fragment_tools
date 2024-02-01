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
> Requires python version **>=3.8** and **<= 3.12**

### From pypi:

{% highlight bash %}
> pip install scatac_fragment_tools
{% endhighlight %}

### From source:

{% highlight bash %}
> git clone https://github.com/aertslab/scatac_fragment_tools.git
> cd scatac_fragment_tools
> pip install.
{% endhighlight %}

## Testing installation

{% highlight bash %}
> scatac_fragment_tools
Usage: scatac_fragment_tools [-h] {bigwig,split} ...

scATAC-fragment-tools (v0.1.0): Tools for processing scATAC-seq fragments.

Options:
  -h, --help      show this help message and exit

Subcommands:
  {bigwig,split}

{% endhighlight %}
