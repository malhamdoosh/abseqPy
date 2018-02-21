# Table of contents
* [Introduction](#introduction)
    * [AbSeq](#abseq)
    * [About](#about)
    * [Wiki](#wiki)
* [Setup](#setup)
    * [Dependencies](#dependencies)
    * [Download](#download)
* [Usage](#usage)
    * [Parameter definitions](#parameter-definitions)
    * [Use cases](#use-cases)
        * [Directory of samples](#directory-of-samples)
        * [Single sample](#single-sample)
            * [Python-only plot](#python-only-plot)
    * [R vs Python plots](#r-vs-python-plots)
* [References](#references)

# Introduction

## AbSeq
**AbSeq** is a quality control pipeline for the construction of antibody libraries, currently running on version 1.1.4

## About
* **AbSeq** is developed by Monther Alhamdoosh and JiaHong Fong
* For comments and suggestions, email m.hamdoosh \<at\> gmail \<dot\> com

## Wiki
<!-- TODO -->

There will be more information on contribution guidelines in the wikipage.

* Writing tests
* Code review
* Other guidelines
* How to run tests

# Setup

## Dependencies
Before proceeding, you should confirm that you have at least these tools readily available:
* [C/C++ compilers](https://gcc.gnu.org/)
* [curl](https://curl.haxx.se/)
* [git](https://git-scm.com/)
* [make](https://en.wikipedia.org/wiki/Make_(software))
* [python](https://www.python.org) v2.7
* [CMake](https://cmake.org/)
* [R](https://cran.r-project.org/) (optional, see below)

By default, AbSeq plots with `Rscript`. If you do not want to install R in your system, you can opt out
by specifying `-rs off`, see [usage](#usage) and [python only plot](#python-only-plot) for more information.

## Download

> Warning: You should be cautious before installing AbSeq! Installing AbSeq also installs a set of python libraries
(and their versions - you can see them in requirements.txt) it uses.
If you do not wish to pollute your python environment, consider using
some kind of virtual environment. It's considered a good practice and is widely used,
see [virtualenv](https://packaging.python.org/guides/installing-using-pip-and-virtualenv/) or
[conda](https://conda.io/docs/user-guide/tasks/manage-environments.html).

> If you already have the datafiles required by IgBlast, now is a good time to set your `$IGDATA` and `$IGBLASTDB`
variables to `/path/to/dir/containing/internal_data_and_optional_data` and `/path/to/igblast/makeblastdb_output` to save
some time installing AbSeq.

To download and install AbSeq:
```bash
$ git clone <insert-abseq-repo-url>
$ cd abseq
$ pip install .
``` 

depending on what your system already has _and your internet speed_, this may take awhile to install.

Folks who care about maintaining (and occasionally using) different versions of AbSeq can create a symlink to `abseq-run` script instead.

Either way, installing AbSeq via `pip install .` will always create an executable in your path called `abseq`.

If you really want to, you can install and 'compile' AbSeq yourself. After cloning this repository (omitting `pip install .`), follow the
instructions in INSTALL.md, then you should create a symlink to `abseq-run`.
  
# Usage
## Parameter definitions
| Parameters                 	| Arguments                                                                                                                                                      	| Help                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  	|
|----------------------------	|----------------------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| `-f1, --file1`               	| /path/to/file_one.fast[aq] or /path/to/directory_of_samples/                                                  	                                                | Path to the first input (FASTQ or FASTA) file. Alternatively, this parameter also accepts a directory of samples. When a directory of samples is provided, all the samples will be analyzed simultaneously.                                                                                                                                                                                                                                                                                        	                                    |
| `-f2, --file2`              	| /path/to/file_two.fastq                                                                                                                                        	| Only specify this if the sample is paired-end. This option is also invalid if  `-f1` was a directory. AbSeq will automatically detect paired-end reads if a directory was provided.                                                                                                                                                                                                                                                                                                                                                   	|
| `-c, --chain`              	| hv, lv, kv [default = hv]                                                                                                                                      	| Chain type. Accepted values are: Heavy chains `hv`, Lambda chains `lv`, and Kappa chains `kv`                                                                                                                                                                                                                                                                                                                                                                                                                                          	|
| `-t, --task`                 	| all, annotate, abundance, diversity, fastqc, productivity, primer, 5utr, rsasimple, rsa, seqlen, secretion, seqlenclass [default = abundance]                  	| The task to perform.<br />`all` conducts `annotate`, `abundance`, `diversity`, `fastqc`, `productivity` and `seqlen`. <br />`annotate` runs IgBLAST to annotate reads. <br />`abundance` annotate and analyze germline abundances. <br />`diversity` analyzes the diversity of clonotypes <br />`fastqc` runs FastQC on the reads <br />`productivity` analyzes the productivity of reads <br />`seqlen` plots the sequence length distributions <br />`seqlenclass` (todo) <br />`primer` (todo) <br />`5utr` (todo) <br />`rsasimple` (todo) <br />`rsa` (todo) <br />`secretion` (todo) 	|
| `-s, --seqtype`              	| dna, protein [default = dna]                                                                                                                                   	| Specify the sequence type.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            	|
| `-m, --merger`               	| flash, leehom [default = leeHom]                                                                                                                               	| Specify the merger used to merge paired-end reads. This option is invalid if  `-f1` is a FASTA file. This merger will also be used for *all* samples if  `-f1` is a directory.                                                                                                                                                                                                                                                                                                                                                        	|
| `-o, -outdir`                	| /path/to/output/directory/ [default = current dir]                                                                                                             	| Output directory for reports and plots.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               	|
| `-n, --name`                 	| analysis_name [default = name of file 1]                                                                                                                       	| Reports and plots will refer to the sample with this name. This option is **ignored** if  `-f1` is a directory.                                                                                                                                                                                                                                                                                                                                                                                                                       	|
| `-cl, --clonelimit`          	| integer value [default = 100]. Use `-cl inf` to retain *all* clones.                                                                                           	| AbSeq reports the overly (and underly) expressed top N clonotypes under the directory `./<outdir>/<name>/diversity/clonotypes/`.                                                                                                                                                                                                                                                                                                                                                                                                      	|
| `-b, --bitscore`             	| range values (inclusive). If only one value is provided, implicitly starts from the value. Example: `-b 300-inf` is equivalent to `-b 300` [default = 0-inf]    	| **Filtering criterion**  <br /> Filters sequences based on the V-gene bitscore from IgBlast.                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `-ss, --sstart`              	| range values (inclusive). If only one value is provided, implicitly starts from the value. Example: `-ss 1-inf` is equivalent to `-ss 1` [default = 1-inf]     	| **Filtering criterion**  <br /> Filters sequences based on the subject's V-gene starting index.                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `-qs, --qstart`              	| range values (inclusive). If only one value is provided, implicitly starts from the value. Example: `-qs 1-inf` is equivalent to `-qs 1` [default = 1-inf]     	| **Filtering criterion**  <br /> Filters sequences based on the query 's V-gene starting index.                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `-al, --alignlen`            	| range values (inclusive). If only one value is provided, implicitly starts from the value. Example: `-al 300-inf` is equivalent to `-al 300` [default = 0-inf] 	| **Filtering criterion**  <br /> Filters sequences based on the sequence alignment length.                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `-qo, --qoffset`             	| integer value [default = automatically inferred for every sequence]                                                                                               |  Subsequences prior to this index are ignored. (Can be used to ignore primer residue in front of sequences). By default, each individual sequence's offset is inferred automatically.                                                                                                                                                                                                                                                                                                                                                     |
| `-u, --upstream`             	| range values. [default = 1-inf]                                                                                                                                	| Range of upstream sequences, secretion signal analysis.                                                                                                                                                                                                                                                                                                                                                                                                                                                                               	|
| `-t5, --trim5`               	| integer value [default = 0]                                                                                                                                    	| Number of nucleotides to trim on the 5'end of V gene                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  	|
| `-t3, --trim3`               	| integer value [default = 0]                                                                                                                                    	| Number of nucleotides to trim on the 3'end of V gene                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  	|
| `-p5off, --primer5endoffset` 	| integer value [default = 0]                                                                                                                                    	| Number of nucleotides for 5'end offset                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                	|
| `-p, --primer`               	| integer value [default = -1]                                                                                                                                   	| Not implemented yet                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   	|
| `-d, --database`             	| /path/to/database/ Defaults to `$IGBLASTDB` env variable setup [earlier](#exporting-environment-variables)                                                     	| `/path/to/igblastDB/databases/` where the output of `makeblastdb` resides                                                                                                                                                                                                                                                                                                                                                                                                                                                               	|
| `-q, --threads`              	| integer value [default = 8]                                                                                                                                    	| Number of processes running concurrently                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              	|
| `-rs, --rscripts`             | empty, string, or a file [default = empty]                                                                                                                        | Reporting engine and sample comparison flag. When `-rs <arg>` is specified, where: <br /> `<arg>` is a string of samples, then there will be R plots generated to compare the samples explicitly. <br /> `<arg>` is a filename, then there will be R plots generated to compare the samples explicitly. The sample comparison are separated by newlines. <br /> `<arg>` is the string `off`, then AbSeq plots in python instead of R. Note that there is no sample comparison available when plotting in python. <br /> `<arg>` is not specified, then the default behaviour of plotting in R with no explicit sample comparison is conducted. <br /> Note that the supplied string or file for sample comparison only makes sense when `-f1` is supplied with a directory. <br /> See [use case examples](#directory-of-samples) for examples of this flag in use.
| `-r, --report-interim`       	| no arguments                                                                                                                                                   	| Specify this flag to generate report. Not implemented yet.                                                                                                                                                                                                                                                                                                                                                                                                                                                                            	|
| `-f4c, -fr4cut`              	| no arguments                                                                                                                                                   	| Specify this flag to remove subsequence after FR4 region. Not specified by default.                                                                                                                                                                                                                                                                                                                                                                                                                                                   	|
| `-st, --sites`               	| /path/to/file                                                                                                                                                  	| Required if `-t rsa` or `-t rsasimple` is specified. (todo)                                                                                                                                                                                                                                                                                                                                                                                                                                                                           	|
| `-p3, --primer3end`          	| /path/to/file                                                                                                                                                  	| Path to primer 3'end file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             	|
| `-p5, --primer5end`          	| /path/to/file                                                                                                                                                  	| Path to primer 5'end file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             	|
| `-v, --version`              	| no arguments                                                                                                                                                   	| Prints the current AbSeq version and exits                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            	|
| `-h, --help`                 	| no arguments                                                                                                                                                   	| Prints this help message and exits                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    	|


## Use cases

### Directory of samples
Let `/Users/john/sequences/` be a directory as such:
```bash
$ ls /Users/john/sequences/
SRR1_BZ123_CAGGG-GGACT_L001.fasta.gz
SRR2_BZ929_CAGGG-GGACT_L001_R1.fastq.gz
SRR2_BZ929_CAGGG-GGACT_L001_R2.fastq.gz 
SRR3_BZ929_CAGGG-GGACT_L001.fastq.gz 
```
We wish to analyze all samples and explicitly compare samples SRR1 and SRR2, and also SRR1, SRR2, and SRR3.

Then running AbSeq:
```bash
$ abseq -f1 /Users/john/sequences/ -t all \
    -o results -m leehom -q 4 \
    -cl inf -ss 1.0-3.0 \
    -b 350 \
    -al 260 \
    -rs "SRR1 | SRR2; SRR1|SRR2_BZ929_L001|SRR3"
```

will produce analysis for:

  * `SRR1_L001` in `results/SRR1_BZ123_CAGGG-GGACT_L001/`,
  * `SRR2_L001` in `results/SRR2_BZ929_CAGGG-GGACT_L001/`,
  * `SRR3_L001` in `results/SRR3_BZ929_CAGGG-GGACT_L001/`,
  
sample comparison analysis for:

  * `SRR1_L001` and `SRR2_L001` in `results/SRR1_L001_vs_SRR2_L001/`,
  * `SRR1_L001`, `SRR2_L001`, and `SRR3_L001` in `results/SRR1_L001_vs_SRR2_L001_vs_SRR3_L001/`.

**Note** that `;` in `-rs` separates different _sample comparison_ while `|` _separates samples_ **within** the same comparison.

> Tip 0: AbSeq only runs analysis on the samples provided in `-rs` (if `-rs` was specified)
rather than _all_ samples in `-f1 dirname`. 

> Tip 1: `-rs` uses fuzzy string search. Therefore, providing either the full sample file name or truncated name will work
fine. Observe in the example, the mixed use of `SRR2_BZ929_L001` and `SRR2`. `SRR2_L001` and `SRR2_BZ929_CAGGG-GGACT_L001_R1.fastq.gz` works just fine too.

> Tip 2: `-rs` ignores whitespaces. Feel free to put spaces between `|`s and `;`s.

Assuming
```bash
$ cat comparison.cfg
SRR1 | SRR2
SRR1_L001 | SRR2 | SRR3
```
> Tip 3: Then `-rs comparison.cfg` would've worked exactly the same as `-rs "SRR1 | SRR2; SRR1|SRR2_BZ929_L001|SRR3"`. This is especially useful when there are long or complicated comparisons!


### Single sample
Running AbSeq on one sample:
```bash
$ abseq -f1 SRR2_BZ929_CAGGG-GGACT_L001_R1.fastq.gz \
        -f2 SRR2_BZ929_CAGGG-GGACT_L001_R2.fastq.gz \
        -t all -o results -m leehom -q 4 \
        -cl inf -ss 1.0-3.0 -b 350 -al 260
```
will produce analysis for `SRR2_L001` in `results/SRR2_BZ929_CAGGG-GGACT_L001/`.

#### Python-only plot
To turn off R plots and use python's plotting engine instead:
```bash
$ abseq -f1 SRR2_BZ929_CAGGG-GGACT_L001_R1.fastq.gz \
        -f2 SRR2_BZ929_CAGGG-GGACT_L001_R2.fastq.gz \
        -t all -o results -m leehom -q 4  \
        -cl inf -ss 1.0-3.0 -b 350  -al 260 \
        -rs off
```

## R vs Python plots
The main difference between R and Python's plots in AbSeq is the ability to compare multiple samples.

If `-f1` is a directory of samples and `-rs` was provided with sample comparisons, the R plots generated by AbSeq includes explicit comparisons between them.

If `-f1` is a directory of samples or a file but `-rs` is supplied with the argument `off`, then there will only be individual sample plots _in python_.

Finally, if `-f1` is a directory and `-rs` is not specified or `-rs` is specified with no aruguments, then the individual sample plots will be 
generated in R but there will be **no** explicit sample comparisons.

# References 

Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD,
Higgins DG (2011). 
[Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega](http://www.nature.com/msb/journal/v7/n1/full/msb201175.html).
Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75


[leeHom](https://github.com/grenaud/leeHom): adaptor trimming and merging for Illumina sequencing reads
Gabriel Renaud, Udo Stenzel, Janet Kelso
Nucleic Acids Res. 2014 Oct 13; 42(18): e141. Published online 2014 Aug 6. doi: 10.1093/nar/gku699
PMCID: PMC4191382

Andrews S. (2010). [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): a quality control tool for
high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

[TAMO](http://fraenkel.mit.edu/TAMO/): a flexible, object-oriented framework for analyzing transcriptional regulation using DNA-sequence motifs. 
Bioinformatics. 2005 Jul 15;21(14):3164-5.

[IgBLAST](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3692102/): an immunoglobulin variable domain sequence analysis tool
Jian Ye, Ning Ma, Thomas L. Madden, James M. Ostell
Nucleic Acids Res. 2013 Jul; 41(Web Server issue): W34–W40. Published online 2013 May 11. doi: 10.1093/nar/gkt382
PMCID: PMC3692102
