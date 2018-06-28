# Mandatory dependecies
* [Clustal Omega](http://www.clustal.org/omega/) v1.2.1 or higher
    - Download and extract the tarball or download and use the pre-compiled binaries
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) v0.11.5 or higher
    - Download and extract the binaries
* [leeHom](https://github.com/grenaud/leeHom)
    - Only one of `leeHom`, [`FLASH`, and `PEAR`](#optional-dependencies) is required.
    - This is the default merger used by AbSeq on non Windows machines. AbSeq switches to [FLASH](#optional-dependencies) as
    the default merger if it detects Windows; Windows users might find it easier to just download the pre-built
    `FLASH` binary.
    - Follow the installation guide in their README , leeHom uses `CMake` and `make` as their build tool.
* [IgBLAST](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/) v1.7
    - Setup guide [here](https://ncbi.github.io/igblast/cook/How-to-set-up.html)
    - Make sure to follow **_every_** step detailed in the guide
    - **_Important_**: The environment variable `$IGDATA` must be exported, as stated in the guide.
* [Ghostscript](https://www.ghostscript.com/download/gsdnld.html) v9.22 or higher
    - Download and follow the instructions to install [here](https://www.ghostscript.com/doc/9.22/Install.htm)

# Optional dependencies
* [FLASH](https://sourceforge.net/projects/flashpage/files/) v1.2.11 or higher
    - Only required if `leeHom` and `PEAR` is not installed
    - Download, extract, `make`, or download the pre-built binary
* [PEAR](https://www.h-its.org/downloads/pear-academic/#release) any version
    - Only required if `FLASH` and `leeHom` is not installed
