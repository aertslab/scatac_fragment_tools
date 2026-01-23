---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
title: split
nav_order: 3
---

# split

Split fragment files by cell type.

## Usage and options

{% highlight bash %}

scatac_fragment_tools split \
    -f <PATH_TO_SAMPLE_TO_FRAGMENT_DEFINITION> \
    -b <PATH_TO_CELL_TYPE_TO_CELL_BARCODE_DEFINITION> \
    -c <CHROM_SIZES_FILENAME> \
    -o <PATH_TO_OUTPUT_FOLDER>

{% endhighlight %}

<b><ins>Required Arguments</ins></b>

**-f, --sample_fragments**
{: .py-0 .text-blue-300}
Path to a text file mapping sample names to fragment files.
{: .px-6 .py-0}

**-b, --cell_type_barcodes**
{: .py-0 .text-blue-300}
Path to a text file mapping samples to cell types and cell types to cell barcodes.
{: .px-6 .py-0}

**-c, --chrom**
{: .py-0 .text-blue-300}
Filename with chromosome sizes (*.chrom.sizes, *.fa.fai).
{: .px-6 .py-0}


**-o, --output**
{: .py-0 .text-blue-300}
Path to output folder.
{: .px-6 .py-0}

<b><ins>Optional arguments</ins></b>

**-t, --temp**
{: .py-0 .text-blue-300}
Path to temporary folder. Default: /tmp
{: .px-6 .py-0}

**-n, --n_cpu**
{: .py-0 .text-blue-300}
Number of cores to use. Default: 1
{: .px-6 .py-0}

**-v, --verbose**
{: .py-0 .text-blue-300}
Whether to print progress. Default: False
{: .px-6 .py-0}

**--clear_temp**
{: .py-0 .text-blue-300}
Whether to clear the temporary folder. Default: False
{: .px-6 .py-0}

**-s, --sep**
{: .py-0 .text-blue-300}
Separator for text files. Default: '\t'
{: .px-6 .py-0}

**--sample_column**
{: .py-0 .text-blue-300}
Column name for the sample name. Default: sample
{: .px-6 .py-0}

**--fragment_column**
{: .py-0 .text-blue-300}
Column name for the path to the fragment file. Default: path_to_fragment_file
{: .px-6 .py-0}

**--cell_type_column**
{: .py-0 .text-blue-300}
Column name for the cell type. Default: cell_type
{: .px-6 .py-0}

**--cell_barcode_column**
{: .py-0 .text-blue-300}
Column name for the cell barcode. Default: cell_barcode
{: .px-6 .py-0}

**--add_sample_id**
{: .py-0 .text-blue-300}
Prefix sample id to cell barcode in pseudobulk fragment file. Default: False
{: .px-6 .py-0}

## Examples of input files

**sample_to_fragment.tsv**

{% highlight bash %}

sample  path_to_fragment_file
A       a.fragments.tsv.gz
B       b.fragments.tsv.gz

{% endhighlight %}

**cell_type_to_cell_barcode.tsv**

{% highlight bash %}

sample  cell_type  cell_barcode
A       type_1     TTAGCTTAGGAGAACA-1
A       type_1     TTAGCTTAGGAGAACA-1
A       type_1     ATATTCCTCTTGTACT-1
A       type_2     TGTGACAGTACAACGG-1
A       type_2     CATGCCTTCTCTGACC-1
A       type_2     ATCGAGTAGGTTCGAG-1
A       type_3     CTCTCAGGTCCCTTTG-1
A       type_3     TTCGGTCTCACGTGTA-1
A       type_3     GTGACATCATTGTTCT-1
A       type_4     AAGGAGCCATCGACCG-1
A       type_4     ACCAAACTCTTAAGCG-1
A       type_4     CATTGGATCTCTTCCT-1
A       type_5     AGGCGAAAGGTCTTTG-1
A       type_5     AACGAGGCATCATGTG-1
A       type_5     CTACTTAGTCATGAGG-1
B       type_1     ATTACCTGTGTGCTTA-1
B       type_1     CATAACGTCGGTTGTA-1
B       type_1     ATGTCTTTCGGTCCGA-1
B       type_2     CAATCCCGTAGCGTTT-1

{% endhighlight %}

