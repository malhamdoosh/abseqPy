# Table of contents
* [Introduction](#introduction)
    * [AbSeq](#abseq)
    * [About](#about)
    * [Wiki](#wiki)
* [Dependencies](#dependencies)
    * [Binary dependencies](#binary-dependencies)
        * [Mandatory dependencies](#mandatory-dependencies)
        * [Optional dependencies](#optional-dependencies)
    * [R library dependencies](#r-library-dependencies)
* [Setup](#setup)
    * [Download](#download)
    * [Installation and configuration](#installation-and-configuration)
    * [Exporting variables](#exporting-variables)
    * [Exporting environment variables](#exporting-environment-variables)
* [Usage](#usage)
    * [Parameter definitions](#parameter-definitions)
    * [Use cases](#use-cases)
        * [Directory of samples](#directory-of-samples)
        * [Single sample](#single-sample)
            * [Python-only plot](#python-only-plot)
    * [R vs Python plots](#r-vs-python-plots)

# Introduction
## AbSeq
**AbSeq** is a quality control pipeline for the construction of antibody libraries, currently running on version 1.1.3

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

# Dependencies
## Binary dependencies
In this section, we assume Unix's `make` build tool is readily available.
If you're (optionally) building for `leeHom` paired-end merger, you will require [`CMake`](https://cmake.org/download/) too.

AbSeq requires a few external packages available in your system, namely:

  * ### Mandatory dependencies

    * [Clustal Omega](http://www.clustal.org/omega/) v1.2.1
        - Download and extract the tarball or install the pre-compiled binaries
        - Follow the installation guide [here](http://www.clustal.org/omega/INSTALL)
    * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.5
        - Download and extract
        - Follow the installation guide [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt)
    * [FLASh](https://sourceforge.net/projects/flashpage/files/) v1.2.11
        - This is (currently) the default merger used by AbSeq. Only one of `FLAsH`, [`leeHom`, and `PEAR`](#optional-dependencies) is required.
        - Download and extract
        - Execute `make` in root directory of FLASh
    * [IgBLAST](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/) v1.6
        - There's an amazing setup guide [here](https://ncbi.github.io/igblast/cook/How-to-set-up.html)
        - Make sure to follow **_every_** step detailed in the guide
        - **_Important_**: Make sure you export the environment variables `$IGBLASTDB` and `$IGDATA`.
         See [here](#exporting-environment-variables)
    * [Ghostscript](https://www.ghostscript.com/download/) v9.22
        - Download and follow the instructions to install [here](https://www.ghostscript.com/doc/9.22/Install.htm)
         
  * ### Optional dependencies
  
    * [TAMO](http://fraenkel.mit.edu/TAMO/) v1.0 is only required if you specify secretion signal analysis or 5'UTR analysis [(`-t secretion` or `-t 5utr`)](#parameter-definitions)
        - Click on "Download the package", extract the tarball
        - Follow the installation guide [here](http://fraenkel.mit.edu/TAMO/INSTALL)
        - When prompted to install databases, you can safely skip them. AbSeq **doesn't** require any of those
    * [leeHom](https://github.com/grenaud/leeHom) any version is only required if `FLASh` and `PEAR` is not installed
        - Follow the installation guide in their README
        - As mentioned earlier, leeHom uses `CMake` and `make` as their build tool.
    * [PEAR](https://www.h-its.org/downloads/pear-academic/#release) any version is only required if `FLASh` and `leeHom` is not installed
        - Follow instructions on their website. __Be sure to read their license agreement before you download their software__.
    
> Make sure the above programs are available in your [`$PATH` variable](#exporting-variables). Pay special
attention to the versions - these versions were used during the development and testing process of AbSeq.


## R library dependencies
The plotting facilities provided by AbSeq uses a few R libraries. You will require `ggplot2`,
`RcolorBrewer`, `circlize`, `reshape2`, `VennDiagram`, and `plyr`.

**These packages will be _automatically installed_ if they can't be located in your system.**

> Alternatively, you can **avoid R and all its dependencies entirely** if you explicitly tell AbSeq to [plot in python only](#python-only-plot). 
See [R vs Python plots](#r-vs-python-plots).

Pat yourself on the back - all the dependencies are installed - you're almost there!

# Setup
Setting up AbSeq is simple as pie:

  1. Download AbSeq
  2. Install
  3. Exporting path variables to locate the binaries you installed [earlier](#binary-dependencies)
  
## Download
`git clone` this repository or manually download this repository.

## Installation and configuration
Before proceeding any further, make sure you have all the [external dependencies](#dependencies)
installed and ready to go. You will require Python v2.7 on your system with [pip](https://pip.pypa.io/en/stable/installing/) installed.

Additionally, you will also need R installed on your machine unless you disable R plotting. See [above](#r-library-dependencies).

AbSeq depends on several python packages, namely:

  1. [Biopython](http://biopython.org/)
  2. [pandas](https://pandas.pydata.org/)
  3. [numpy](http://www.numpy.org/)
  4. [Weblogo](https://github.com/WebLogo/weblogo) 

To automatically install the above packages, run:
```bash
$ pip install -r requirements.txt
```
in AbSeq's root directory (where requirements.txt lives).

## Exporting variables
To make the [installed binaries](#binary-dependencies) available in your `$PATH` variable:
```bash
export PATH="/path/to/fastqc:/path/to/leehom/:$PATH"
```
in your `.bashrc` or equivalent.

`/path/to/binaries` is the absolute path to the directory where your programs (listed above) are installed,
each separated by colons. Repeat this for every dependency. Alternatively, move all binaries into one
folder (eg, `/Users/john/bin/`) and `export PATH="/Users/john/bin/:$PATH"`

## Exporting environment variables
Make sure `$IGBLASTDB` and `$IGDATA` are exported.
```bash
export IGBLASTDB="/path/to/data/igblastDB/databases/"
export IGDATA="/path/to/data/igblastDB/data/"
```
Again, `/path/to/data/` is the absolute path to the directory where `/igblastDB/ ... /` lives.



# Usage
## Parameter definitions
| Parameters                 	| Arguments                                                                                                                                                      	| Help                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  	|
|----------------------------	|----------------------------------------------------------------------------------------------------------------------------------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| `-f1, --file1`               	| /path/to/file_one.fast[aq] or /path/to/directory_of_samples/                                                  	                                                | Path to the first input (FASTQ or FASTA) file. Alternatively, this parameter also accepts a directory of samples. When a directory of samples is provided, all the samples will be analyzed simultaneously.                                                                                                                                                                                                                                                                                        	                                    |
| `-f2, --file2`              	| /path/to/file_two.fastq                                                                                                                                        	| Only specify this if the sample is paired-end. This option is also invalid if  `-f1` was a directory. AbSeq will automatically detect paired-end reads if a directory was provided.                                                                                                                                                                                                                                                                                                                                                   	|
| `-c, --chain`              	| hv, lv, kv [default = hv]                                                                                                                                      	| Chain type. Accepted values are: Heavy chains `hv`, Lambda chains `lv`, and Kappa chains `kv`                                                                                                                                                                                                                                                                                                                                                                                                                                          	|
| `-t, --task`                 	| all, annotate, abundance, diversity, fastqc, productivity, primer, 5utr, rsasimple, rsa, seqlen, secretion, seqlenclass [default = abundance]                  	| The task to perform.<br />`all` conducts `annotate`, `abundance`, `diversity`, `fastqc`, `productivity` and `seqlen`. <br />`annotate` runs IgBLAST to annotate reads. <br />`abundance` annotate and analyze germline abundances. <br />`diversity` analyzes the diversity of clonotypes <br />`fastqc` runs FastQC on the reads <br />`productivity` analyzes the productivity of reads <br />`seqlen` plots the sequence length distributions <br />`seqlenclass` (todo) <br />`primer` (todo) <br />`5utr` (todo) <br />`rsasimple` (todo) <br />`rsa` (todo) <br />`secretion` (todo) 	|
| `-s, --seqtype`              	| dna, protein [default = dna]                                                                                                                                   	| Specify the sequence type.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            	|
| `-m, --merger`               	| flash, leehom [default = flash]                                                                                                                                	| Specify the merger used to merge paired-end reads. This option is invalid if  `-f1` is a FASTA file. This merger will also be used for *all* samples if  `-f1` is a directory.                                                                                                                                                                                                                                                                                                                                                        	|
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
| `-d, --database`             	| /path/to/database/ Defaults to `$IGBLASTDB` env variable setup [earlier](#exporting-environment-variables)                                                     	| Path to `igblastDB/databases/`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        	|
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
