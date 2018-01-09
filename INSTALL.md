# Installation

> This section provides a comprehensive installation guide, including downloading and installing AbSeq
and all its dependencies without root access. 

To reiterate:
#### Mandatory dependencies
* [Clustal Omega](http://www.clustal.org/omega/) v1.2.1
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.5
* [FLASh](https://sourceforge.net/projects/flashpage/files/) v1.2.11
* [IgBLAST](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/) v1.6
* [Ghostscript](https://www.ghostscript.com/download/) v9.22
         
#### Optional dependencies
* [TAMO](http://fraenkel.mit.edu/TAMO/) v1.0 is only required if you specify secretion signal analysis or 5'UTR analysis [(`-t secretion` or `-t 5utr`)](#parameter-definitions)
* [leeHom](https://github.com/grenaud/leeHom) any version is only required if `FLASh` and `PEAR` is not installed
* [PEAR](https://www.h-its.org/downloads/pear-academic/#release) any version is only required if `FLASh` and `leeHom` is not installed


# 1. Prepare a directory with read/write access
If you don't mind installing binaries into protected areas such as `/usr/local/` and the likes,
you can safely skip this step. You **must** also ignore typing anything related to `Users/john/mirror` whenever
you see it in the following steps. That is, for example:
`bash$ ./configure --prefix=/Users/john/mirror` should be
`bash$ ./configure` instead.

Unfortunate folks with no root access:

For the remainder of this setup, we assume a directory with full read and write access at:
```bash
bash$ pwd
/Users/john/
```

We can then prepare a directory for most of the installation process:
```bash
bash$ mkdir -p mirror/bin
bash$ ls -F
mirror/
```
I called it mirror because it "mirrors" `/usr/local` as we shall see later.

Now, we need to access binaries from `/Users/john/mirror/bin`, we can do this by adding
the following line in your `.bashrc` or `.bash_profile`:
```bash
export PATH="$PATH:/Users/john/mirror/bin/"
```

# 2. Install build tools.
You will require a C/C++ compiler. Get it from your appropriate vendor.

Additionally, you'll require `make` and `CMake`. Typically `make` is shipped with
Unix-like distributions.

Download the compiled [CMake](https://cmake.org/download) for your operating system and copy
the cmake binary over to `/Users/john/mirror/bin`

> CMake is only needed if you want to compile [leeHom](https://github.com/grenaud/leeHom) from source

# 3. Build, compile sources, and install python packages
* Download [argtable2](http://argtable.sourceforge.net/) before [Clustal Omega](http://www.clustal.org/omega/)

`Clustal Omega` depends on `argtable2`, so we shall start with `argtable2`; unpack and `cd` into
`argtable2/`:
```bash
bash$ ./configure --prefix=/Users/john/mirror
bash$ make
bash$ make install
```
Now we can compile `Clustal Omega`; after unpacking and 'cd-ing' into `clustal-omega/`
```bash
bash$ ./configure --prefix=/Users/john/mirror CFLAGS='-I/Users/john/mirror/include' LDFLAGS='-L/Users/john/mirror/lib'
bash$ make
bash$ make install
```
> users installing as root can just type `./configure`

* Download [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), unpack, and copy the `fastqc`
binary over to `/Users/john/mirror/bin`.

Remember to flip the executable bit on: `chmod +x fastqc`

> users installing as root can create a sym link to /usr/local/bin or edit the PATH variable, or copying fastqc to
/usr/local/bin.

* Download [FLASh](https://sourceforge.net/projects/flashpage/files/), unpack, and `cd` into the
directory. Simply type `make` and copy the `flash` binary over to `/Users/john/mirror/bin`

* Download [IgBLAST](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/). They offer
pre-compiled binaries for your operating system. Next, follow the pre-processing steps for
`IgBlast`'s internal data [here](https://ncbi.github.io/igblast/cook/How-to-set-up.html).

* Download [Ghostscript](https://www.ghostscript.com/download/), unpack, and `cd` - you know the drill.
```bash
bash$ ./configure --prefix=/Users/john/mirror/
bash$ make
bash$ make install
```

> Ghostscript is required by WebLogo (weblogo is installed in requirements.txt)

* Install python packages

In AbSeq's root directory (where requirements.txt is):
```bash
bash$ pip install -r requirements.txt
```

