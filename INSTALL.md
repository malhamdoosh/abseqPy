# Installation

> This section provides a comprehensive installation guide, including downloading and installing AbSeq
and all its dependencies without root access. 

If you're here because `pip` wants to install AbSeq with root access, have a look [here](https://stackoverflow.com/a/7465532).

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
    * [leeHom](https://github.com/grenaud/leeHom) any version
        - This is (currently) the default merger used by AbSeq. Only one of `leeHom`, [`FLASh`, and `PEAR`](#optional-dependencies) is required.
        - Follow the installation guide in their README
        - As mentioned earlier, leeHom uses `CMake` and `make` as their build tool.
    * [IgBLAST](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/) v1.8
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
    * [FLASh](https://sourceforge.net/projects/flashpage/files/) v1.2.11 only required if `leeHom` and `PEAR` is not installed
        - Download and extract
        - Execute `make` in root directory of FLASh
    * [PEAR](https://www.h-its.org/downloads/pear-academic/#release) any version is only required if `FLASh` and `leeHom` is not installed
        - Follow instructions on their website. __Be sure to read their license agreement before you download their software__.
    
> Make sure the above programs are available in your [`$PATH` variable](#exporting-variables). Pay special
attention to the versions - these versions were used during the development and testing process of AbSeq.


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

`Clustal Omega` depends on `argtable2`, so we shall start with `argtable2`; unpack and `cd` into `argtable2/`:
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

## R library dependencies
The plotting facilities provided by AbSeq uses a few R libraries. You will require `ggplot2`,
`RcolorBrewer`, `circlize`, `reshape2`, `VennDiagram`, and `plyr`.

**These packages will be _automatically installed_ if they can't be located in your system.**

> Alternatively, you can **avoid R and all its dependencies entirely** if you explicitly tell AbSeq to plot in python only.
See R vs Python plots in README.

## Exporting variables
If you chose to install the binaries in some other location, you should make them visible in your `$PATH` variable.
(You can ignore this if you've followed the instructions above to a tee, we've already done this when
we appended mirror/bin to our `PATH`)

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
where `igblastDB/databases/` should have all the processed files required for `IgBlast` (this is the output of `makeblastdb`)
and `igblastDB/data/` should contain 2 sub-folders (you can obtain them from IgBlast's ftp site).
```bash
$ ls /path/to/data/igblastDB/databases/
imgt_human_ighc          imgt_human_igkc          imgt_human_igkv_p.pin    imgt_human_iglv.nsd      imgt_mouse_ighv.nsd
imgt_human_ighc_p        imgt_human_igkc.nhr      imgt_human_igkv_p.pog    imgt_human_iglv.nsi      imgt_mouse_ighv.nsi
imgt_human_ighd          imgt_human_igkc.nin      imgt_human_igkv_p.psd    imgt_human_iglv.nsq      imgt_mouse_ighv.nsq
imgt_human_ighd.nhr      imgt_human_igkc.nog      imgt_human_igkv_p.psi    imgt_human_iglv_p        imgt_mouse_ighv_p
imgt_human_ighd.nin      imgt_human_igkc.nsd      imgt_human_igkv_p.psq    imgt_human_iglv_p.phr    imgt_mouse_ighv_p.phr
imgt_human_ighd.nog      imgt_human_igkc.nsi      imgt_human_iglc          imgt_human_iglv_p.pin    imgt_mouse_ighv_p.pin
imgt_human_ighd.nsd      imgt_human_igkc.nsq  < more entries omited>

$ ls /path/to/data/igblastDB/data/
internal_data/ optional_file/
```

