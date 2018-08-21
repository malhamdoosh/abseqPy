# abseqPy's external dependencies

This document describes the required dependencies to get `abseqPy` working. There is a handy
[installation script](install.py) that automates this process, but this document serves as a
backup if manual installation is required.

## Mandatory dependencies 
* [Clustal Omega](http://www.clustal.org/omega/) v1.2.1 or higher
    - Download and extract the tarball or download and use the pre-compiled binaries
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) v0.11.5 or higher
    - Download and extract the binaries
* [leeHom](https://github.com/grenaud/leeHom) latest GitHub version
    - Only one of `leeHom`, [`FLASH`, `PEAR`](#optional-dependencies) is required.
    - This is the default merger used by AbSeq on non Windows machines. AbSeq switches to [FLASH](#optional-dependencies) as
    the default merger if it detects Windows; Windows users might find it easier to just download the pre-built
    `FLASH` binary.
    - Follow the installation guide in their README , leeHom uses `CMake` and `make` as their build tool.
* [IgBLAST](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/) v1.7
    - Make sure to follow **_every_** step detailed in the [guide](https://ncbi.github.io/igblast/cook/How-to-set-up.html)
    - **_Important_**: The environment variable `$IGDATA` must be exported, as stated in the guide.
* [Ghostscript](https://www.ghostscript.com/download/gsdnld.html) v9.22 or higher
    - Download and follow the instructions to install [here](https://www.ghostscript.com/doc/9.22/Install.htm)

## Optional dependencies
* [FLASH](https://sourceforge.net/projects/flashpage/files/) v1.2.11 or higher
    - Required if `leeHom` and `PEAR` is not installed
    - Download, extract, `make` or download the pre-built binary
* [PEAR](https://www.h-its.org/downloads/pear-academic/#release) any version
    - Required if `FLASH` and `leeHom` is not installed
    

## IgBLAST configuration

A final note, if the [installation script](install.py) was not used, be sure to set the environment variable `IGBLASTDB` to
the path where the germline V, D, and J gene sequences are. (i.e. the directory where `my_seq_file` lives in the
example in IgBLAST's [how to setup](https://ncbi.github.io/igblast/cook/How-to-set-up.html))

If `IGBLASTDB` is not set, `abseqPy` will require the `-d` or `--database` flag (see `abseq -h`)
to be specified with the same directory as `IGBLASTDB` would have had.
