---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
title: bigwig
nav_order: 2
---

# bigwig

Calculate genome coverage for fragments and write result to bigWig file.

## Usage and options

{% highlight bash %}

scatac_fragment_tools bigwig -i <FRAGMENTS_FILENAME> -c <CHROM_SIZES_FILENAME> -o <BIGWIG_FILENAME>

{% endhighlight %}

<b><ins>Required arguments</ins></b>

**-c, --chrom**
{: .py-0 .text-blue-300}
Filename with chromosome sizes (*.chrom.sizes, *.fa.fai).
{: .px-6 .py-0}

**-i, --frag**
{: .py-0 .text-blue-300}
Fragments TSV file for which to calculate genome coverage.
{: .px-6 .py-0}

**-o, --bw**
{: .py-0 .text-blue-300}
BigWig filename to which the genome coverage data will be written.
{: .px-6 .py-0}

<b><ins>Optional arguments</ins></b>

**-n, --normalize**
{: .py-0 .text-blue-300}
Normalize genome coverage data by dividing by sequencing depth (number of fragments) multiplied by 1 million. Default: False
{: .px-6 .py-0}

**-s, --scaling**
{: .py-0 .text-blue-300}
Scaling factor for genome coverage data. If normalization is enabled, scaling is applied afterwards. Default: 1.0
{: .px-6 .py-0}

**-x, --cut-sites**
{: .py-0 .text-blue-300}
Use 1 bp Tn5 cut sites (start and end of each fragment) instead of whole fragment length for coverage calculation. Default: False
{: .px-6 .py-0}

**--chrom-prefix**
{: .py-0 .text-blue-300}
Add chromosome prefix to each chromosome name found in the fragments file. Default: False
{: .px-6 .py-0}
