# Table of contents
* [Introduction](#introduction)
    * [AbSeq](#abseq)
    * [About](#about)
    * [Wiki](#wiki)
    * [Supported platforms](#supported-platforms)
* [Setup](#setup)
    * [Dependencies](#dependencies)
    * [AbSeqPy](#abseqpy)
* [Usage](#usage)
    * [Help](#help)
    * [YAML](#yaml)

# Introduction

## AbSeq
**AbSeq** is a quality control pipeline for high throughput antibody sequencing, currently running on version 1.0

## About
* **AbSeq** is developed by Monther Alhamdoosh and JiaHong Fong
* For comments and suggestions, email m.hamdoosh \<at\> gmail \<dot\> com

## Supported platforms

**AbSeq** has been tested on Unix-like platforms (`macOS` and `Linux`).
Windows support hasn't been tested but should work as expected.

Currently, **AbSeq** runs on Python 2.7 &mdash; fret not, Python 3.6 support is underway.


## Wiki
<!-- TODO -->

There will be more information on contribution guidelines in the wikipage in the future.

* Writing tests
* Code review
* Other guidelines
* How to run tests

# Setup

**AbSeq** depends on external software to work and they should be properly installed and configured before running **AbSeq**.

## Dependencies

A python script is available [here](install.py) which installs all the necessary external dependencies.

Before proceeding, these tools must be present:

* [C/C++ compilers](https://gcc.gnu.org/) (Windows users will not require this)
* [perl](https://www.perl.org/get.html)
* [git](https://git-scm.com/)
* [make](https://en.wikipedia.org/wiki/Make_(software))
* [python](https://www.python.org)
* [CMake](https://cmake.org/)

then, installing external dependencies into a folder named `~/.local/abseq` is as easy as

```bash
$ python install.py ~/.local/abseq
```

The script works with either Python versions, and `~/.local/abseq` can be replaced with any directory.

> The script will remind users to update their environment variables

If the script fails, the list of software required and their version are also listed [here](extdeps.md).



## AbSeqPy

```bash
$ git clone <abseq url>
$ cd abseqPy
$ pip install .
```

`abseq` should now be available on your command line.

> **AbSeq** also installs other python packages, consider using a python virtual environment to prevent overriding existing packages. See [virtualenv](https://packaging.python.org/guides/installing-using-pip-and-virtualenv/) or
[conda](https://conda.io/docs/user-guide/tasks/manage-environments.html).


# Usage
## Help

Invoking `abseq -h` in the command line will display the arguments **AbSeq** uses.


## YAML

Besides calling `abseq` with command line arguments, `abseq` also supports `-y <file>` or `--yaml <file>` argument that
reads off arguments defined within the provided `file`.

#### Example
Assuming a file named `example1.yml` has the following content:

```yaml
defaults:
    task: all
    outdir: results
    threads: 7
    bitscore: 350
    alignlen: 260
    sstart: 1-3
---
name: PCR1
file1: fastq/PCR1_ACGT_L001_R1.fastq.gz
file2: fastq/PCR1_ACGT_L001_R2.fastq.gz
---
name: PCR2
file1: fastq/PCR2_ACGT_L001_R1.fastq.gz
file2: fastq/PCR2_ACGT_L001_R2.fastq.gz
---
name: PCR3
file1: fastq/PCR3_ACGT_L001_R1.fastq.gz
file2: fastq/PCR3_ACGT_L001_R2.fastq.gz
bitscore: 300
task: abundance
detailedComposition: ~
```

then executing `abseq -y example1.yml` is equivalent to running 3 simultaneous instances of
`abseq` with the arguments in the `defaults` field applied to each sample (PCR1, PCR2, PCR3):

```bash
$ abseq --task all --outdir results --threads 7 --bitscore 350 --alignlen 260 --sstart 1-3 \
>   --name PCR1 --file1 fastq/PCR1_ACGT_L001_R1.fastq.gz --file2 fastq/PCR1_ACGT_L001_R2.fastq.gz
$ abseq --task all --outdir results --threads 7 --bitscore 350 --alignlen 260 --sstart 1-3 \
>   --name PCR2 --file1 fastq/PCR2_ACGT_L001_R1.fastq.gz --file2 fastq/PCR2_ACGT_L001_R2.fastq.gz
$ abseq --task abundance --outdir results --threads 7 --bitscore 300 --alignlen 260 --sstart 1-3 \
>   --name PCR3 --detailedComposition --file1 fastq/PCR3_ACGT_L001_R1.fastq.gz --file2 fastq/PCR3_ACGT_L001_R2.fastq.gz
```

> The sample *PCR3* overrides the `bitscore` and `task` argument with `300` and `abundance`, and enables the `--detailedComposition` flag
