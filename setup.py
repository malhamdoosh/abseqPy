from __future__ import print_function

import os
import sys
import pip
import shutil
import glob
import tarfile
import platform
import zipfile
import struct
import subprocess
import re

from subprocess import check_output
from ftplib import FTP
from distutils.version import LooseVersion
from setuptools import setup, find_packages
from setuptools.command.install import install

MAC = 'Darwin'
LIN = 'Linux'
WIN = 'Windows'

# VERSIONING:
# 1. [singleton] ==> minimum version
# 2. [min, max]  ==> accepted range of versions
# 3. True        ==> any version
versions = {
    'clustalo': ["1.2.2", '1.2.4'],
    'leehom': True,
    'flash': True,
    'igblast': ['1.7.0'],
    'fastqc': ['0.11.6', '0.11.7'],
    'gs': ['9.22']
}


class FTPBlast:
    def __init__(self, addr, version):
        self.ftp = FTP(addr)
        self.ftp.login()
        self.version = version

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.ftp.quit()

    def install_bins(self, binary, installation_dir):
        path = '/blast/executables/igblast/release/{}/'.format(self.version) + binary
        installation_path = os.path.join(installation_dir, binary)
        if not os.path.exists(installation_dir):
            os.makedirs(installation_dir)
        with open(installation_path, "wb") as fp:
            self.ftp.retrbinary('RETR ' + path, fp.write)

        old_dir = os.path.abspath(".")
        os.chdir(installation_dir)
        tar = tarfile.open(binary, "r:gz")
        tar.extractall()
        tar.close()
        os.chdir(old_dir)

        return glob.glob(os.path.join(installation_dir, 'ncbi-igblast-' + self.version, 'bin') + os.path.sep + '*')

    def download_edit_imgt_pl(self, download_dir):
        path = '/blast/executables/igblast/release/edit_imgt_file.pl'
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)
        download_path = os.path.join(download_dir, 'edit_imgt_file.pl')
        with open(download_path, "wb") as fp:
            self.ftp.retrbinary('RETR ' + path, fp.write)
        os.chmod(download_path, 0o777)

    def download_internal_data(self, download_dir):
        path = '/blast/executables/igblast/release/internal_data/'
        self.ftp.cwd(path)
        species = self.ftp.nlst()
        for s in species:
            self.ftp.cwd(s)
            filenames = self.ftp.nlst()
            download_path = os.path.join(download_dir, 'internal_data', s)
            os.makedirs(download_path)
            for filename in filenames:
                # ignore rhesus_monkey's CVS directory
                if filename == 'CVS':
                    continue
                with open(os.path.join(download_path, filename), "wb") as fp:
                    self.ftp.retrbinary('RETR ' + filename, fp.write)
            self.ftp.cwd('../')

    def download_optional_file(self, download_dir):
        path = '/blast/executables/igblast/release/optional_file/'
        self.ftp.cwd(path)
        filenames = self.ftp.nlst()
        download_path = os.path.join(download_dir, 'optional_file')
        os.makedirs(download_path)
        for filename in filenames:
            with open(os.path.join(download_path, filename), "wb") as fp:
                self.ftp.retrbinary('RETR ' + filename, fp.write)


def _get_sys_info():
    return platform.system(), 8 * struct.calcsize("P")


def _get_software_version(prog):
    try:
        if prog == 'igblast':
            retval = check_output(['igblastn', '-version']).split('\n')[1].strip().split()[2].rstrip(',')
            return retval
        elif prog == 'clustalo' or prog == 'fastqc' or prog == 'gs':
            retval = check_output([prog, '--version']).strip()
            if prog == 'fastqc':
                retval = retval.split()[-1].strip().lstrip("v")
            return retval
        elif prog == 'leehom':
            # leehomMulti, any version
            check_output(['which', 'leeHomMulti'])
        elif prog == 'flash':
            # flash, any version
            check_output(['which', 'flash'])
    except (subprocess.CalledProcessError, OSError):
        return False
    return True


def _needs_installation(prog):
    v = versions[prog]
    software_version = _get_software_version(prog)
    if type(v) == bool or type(software_version) == bool:
        return software_version != v
    if type(v) == list:
        if len(v) == 1:
            return LooseVersion(software_version) < LooseVersion(v[0])
        elif len(v) == 2:
            return not (LooseVersion(v[0]) <= LooseVersion(software_version) <= LooseVersion(v[1]))
        else:
            _error("Unknown versioning scheme")


def _error(msg, stream=sys.stderr, abort=1):
    print(msg, file=stream)
    if abort:
        sys.exit(abort)


def _syml(src, dest):
    plat, _ = _get_sys_info()
    if plat == WIN:
        return
    if not os.path.exists(dest):
        os.makedirs(dest)
    binary_name = os.path.basename(src)
    if src:
        link_src = os.path.abspath(src)
        link_dest = os.path.join(dest, binary_name)
        # anaconda / conda doesn't like os.symlink
        if 'continuum' in sys.version.lower() or 'anaconda' in sys.version.lower():
            _ = check_output(['ln', '-s', link_src, link_dest])
        else:
            os.symlink(link_src, link_dest)


def setup_dir(root):
    from abseq.config import EXTERNAL_DEP_DIR
    output = os.path.join(root, EXTERNAL_DEP_DIR)
    if os.path.exists(output):
        _error("{} already exists! Remove the directory and try again".format(EXTERNAL_DEP_DIR))
    return output


def install_clustal_omega(installation_dir=".", version=versions['clustalo'][-1]):
    # can't use versions yet, pre-compiled binaries are a little out of sync

    from six.moves.urllib import request
    plat, bit = _get_sys_info()
    # clustalo needs to create a dir
    clustalo_installation_dir = os.path.join(installation_dir, 'clustal-omega')
    if not os.path.exists(clustalo_installation_dir):
        os.makedirs(clustalo_installation_dir)

    binary = os.path.join(clustalo_installation_dir, 'clustalo')
    if plat == MAC:
        addr = 'http://www.clustal.org/omega/clustal-omega-1.2.3-macosx'
        request.urlretrieve(addr, binary)
        # add execution bit
        os.chmod(binary, 0o777)
    elif plat == LIN:
        if bit == 64:
            addr = 'http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64'
        elif bit == 32:
            addr = 'http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-32-bit'
        else:
            _error('Unknown architecture. Detected a non 32 or 64 bit system.')
        # noinspection PyUnboundLocalVariable
        request.urlretrieve(addr, binary)
        # add execution bit
        os.chmod(binary, 0o777)
    elif plat == WIN:
        windows_bin = 'clustal-omega-1.2.2-win64.zip'
        addr = 'http://www.clustal.org/omega/' + windows_bin
        request.urlretrieve(addr, windows_bin)
        zip_ref = zipfile.ZipFile(windows_bin)
        zip_ref.extractall(clustalo_installation_dir)
        zip_ref.close()
        # windows has no symlink, go straight to bin directory!
        for f in os.listdir(os.path.join(clustalo_installation_dir, windows_bin[:windows_bin.find('.zip')])):
            src_ = os.path.join(clustalo_installation_dir, windows_bin[:windows_bin.find('.zip')], f)
            shutil.move(src_, os.path.join(installation_dir, 'bin'))
    else:
        _error('Unknown system architecture. Non windows, mac or linux detected')

    return binary


def install_fastqc(installation_dir=".", version=versions['fastqc'][-1]):
    from six.moves.urllib import request
    plat, _ = _get_sys_info()
    addr = 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v{}.zip'.format(version)
    zipname = os.path.join(installation_dir, os.path.basename(addr).strip())
    request.urlretrieve(addr, zipname)
    unzipped_name = 'FastQC'
    zip_ref = zipfile.ZipFile(zipname, 'r')
    zip_ref.extractall(installation_dir)
    zip_ref.close()
    fastqc_dir = os.path.join(installation_dir, unzipped_name)
    # windows has no symlinks, move to bin immediately
    if plat == WIN:
        for f in os.listdir(fastqc_dir):
            shutil.move(os.path.join(fastqc_dir, f), os.path.join(installation_dir, 'bin', f))
        binary = os.path.join(installation_dir, 'bin', 'fastqc')
    else:
        binary = os.path.join(fastqc_dir, 'fastqc')
        os.chmod(binary, 0o777)
    return binary


def install_leehom(installation_dir='.'):
    addr = 'https://github.com/grenaud/leeHom.git'
    old_dir = os.path.abspath(".")

    # clone into installation dir
    os.chdir(installation_dir)

    # clone
    _ = check_output(['git', 'clone', '--recursive', addr])

    # repo is under 'leeHom'
    os.chdir('leeHom')

    # DO NOT USE -j N here, it's not optimized to use concurrent compilation
    _ = check_output(['make'])

    # go back to our original dir
    os.chdir(old_dir)

    return os.path.join(installation_dir, 'leeHom', 'src', 'leeHomMulti')


def install_flash(installation_dir='.'):
    from six.moves.urllib import request
    addr = "http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11-windows-bin.zip"
    flash_ins_dir = os.path.join(installation_dir, "flash")
    if not os.path.exists(flash_ins_dir):
        os.makedirs(flash_ins_dir)
    flash_zip = 'flash.zip'
    request.urlretrieve(addr, flash_zip)
    zip_ref = zipfile.ZipFile(flash_zip)
    zip_ref.extractall(flash_ins_dir)
    for f in os.listdir(flash_ins_dir):
        shutil.move(os.path.join(flash_ins_dir, f), os.path.join(installation_dir, 'bin', f))


def install_ghost_script(installation_dir='.', threads=2, version=versions['gs'][-1]):
    from six.moves.urllib import request
    plat, bit = _get_sys_info()
    target_dir = os.path.abspath(installation_dir)

    if plat != WIN:
        addr = 'https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs{}/ghostpdl-{}.tar.gz'.format(
            version.replace('.', ''), version)
        tarname = os.path.basename(addr)
        old_dir = os.path.abspath('.')

        os.chdir(installation_dir)
        request.urlretrieve(addr, tarname)
        _ = check_output(['tar', '-xvzf', tarname])
        ghs_dir = os.path.splitext(os.path.splitext(tarname)[0])[0]
        os.chdir(ghs_dir)
        _ = check_output(['./configure', '--prefix={}'.format(target_dir)])
        _ = check_output(['make', '-j', str(threads)])
        _ = check_output(['make', 'install'])
        os.chdir(old_dir)
    else:
        binary = "gs{}w64.exe".format(version.replace('.', ''))
        addr = "http://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs{0}/gs{0}w64.exe"\
            .format(version.replace('.', ''))
        request.urlretrieve(addr, binary)
        os.rename(binary, os.path.join(target_dir, 'bin', 'gs'))
    # dont need to return binary directory, it's already in installation_dir/bin


def install_igblast(installation_dir='.', version=versions['igblast'][-1]):
    plat, _ = _get_sys_info()

    with FTPBlast('ftp.ncbi.nih.gov', version) as blast:
        if plat == MAC:
            bins = blast.install_bins('ncbi-igblast-{}-x64-macosx.tar.gz'.format(version), installation_dir)
        elif plat == WIN:
            bins = blast.install_bins('ncbi-igblast-{}-x64-win64.tar.gz'.format(version), installation_dir)
            for bin_ in bins:
                shutil.move(bin_, os.path.join(installation_dir, 'bin'))
        elif plat == LIN:
            bins = blast.install_bins('ncbi-igblast-{}-x64-linux.tar.gz'.format(version), installation_dir)
        else:
            _error("Unknown platform detected")
    # noinspection PyUnboundLocalVariable
    return bins


def install_TAMO():
    # TAMO comes packed with AbSeq, just need to install it!
    tar = tarfile.open('TAMO.tar.gz', "r:gz")
    tar.extractall()
    tar.close()
    old_dir = os.path.abspath(".")
    os.chdir("TAMO-1.0_120321")
    # install!
    _ = check_output(['python', 'setup.py', 'install'])
    # remove files (tar?)
    os.chdir(old_dir)


def download_imgt(download_dir, species, species_layman):
    from six.moves.urllib import request
    links = [
        "http://www.imgt.org/genedb/GENElect?query=7.14+IGHV&species=",
        "http://www.imgt.org/genedb/GENElect?query=7.14+IGHD&species=",
        "http://www.imgt.org/genedb/GENElect?query=7.14+IGHJ&species=",
        "http://www.imgt.org/genedb/GENElect?query=7.14+IGKV&species=",
        "http://www.imgt.org/genedb/GENElect?query=7.14+IGKJ&species=",
        "http://www.imgt.org/genedb/GENElect?query=7.14+IGLV&species=",
        "http://www.imgt.org/genedb/GENElect?query=7.14+IGLJ&species="
    ]

    path = os.path.join(download_dir, 'imgt_' + species_layman)
    os.makedirs(path)
    for url in links:
        gene = url[url.find("+") + 1:url.find("&")].lower()
        output = "{}_{}.imgt.raw".format(os.path.join(path, species_layman), gene)
        request.urlretrieve(url+species, output)
        with open(output[:output.rfind(".")], "w") as writer, \
                open(output) as reader:

            line = reader.readline()
            while not line.startswith("<b>Number of results"):
                line = reader.readline()
            if not line:
                raise Exception("File has no IMGT sequences")

            reader.readline()  # \n
            reader.readline()  # <pre>

            for line in reader:
                if line.startswith("</pre>"):
                    # finish writing sequences
                    break
                writer.write(line)
            # remove raw
            os.remove(output)


def igblast_compat(edit_imgt_bin, make_blast_bin, data_dir, output_dir):
    from Bio import SeqIO
    for f in os.listdir(data_dir):
        clean_fasta = os.path.join(output_dir, 'imgt_' + f[:f.find(".")])
        with open(clean_fasta, 'w') as fp:
            ret = subprocess.call(['perl', edit_imgt_bin, os.path.join(data_dir, f)], stdout=fp)
            assert ret == 0
        records = []
        seen = {}
        for rec in SeqIO.parse(clean_fasta, 'fasta'):
            rec.description = ''
            rec.seq = rec.seq.upper()
            if rec.id not in seen:
                seen[rec.id] = 0
            else:
                # IGHV1-45*03, IGHV1-45*03_1, IGHV1-45*03_2 ...
                seen[rec.id] += 1
                rec.id = "{}_{}".format(rec.id, str(seen[rec.id]))
            records.append(rec)
                
        SeqIO.write(records, clean_fasta, 'fasta')
        _ = check_output([make_blast_bin, '-parse_seqids', '-dbtype', 'nucl', '-in', clean_fasta])
        if len(re.findall('ig[hkl][vc]', clean_fasta)) > 0:
            clean_fasta_prot = clean_fasta + "_p"
            records = []
            wrong = 0
            for rec in SeqIO.parse(clean_fasta, 'fasta'):
                rec.description = ''
                prot = rec.seq.translate()
                if '*' in str(prot):
                    if str(prot).index('*') == len(prot) - 1:
                        prot = prot[:-1]
                    else:
                        wrong += 1
                        continue
                rec.seq = prot
                records.append(rec)
            SeqIO.write(records, clean_fasta_prot, 'fasta')
            if wrong > 0:
                print('Number of invalid V genes in ' + f + ' is ' + str(wrong) + ' (ignored in the protein db)')
            _ = check_output([make_blast_bin, '-parse_seqids', '-dbtype', 'prot', '-in', clean_fasta_prot])
        print(f + ' has been processed.')


class ExternalDependencyInstaller(install):
    def run(self):
        # although setup() has this, it's installed locally in abseq's installation dir.
        # By pip.installing here, it's going to be available globally
        setup_requires = ['numpy>=1.11.3', 'pytz', 'python-dateutil', 'psutil', 'biopython>=1.66', 'six']
        for pack in setup_requires:
            pip.main(['install', pack])
        # create external deps dir
        d = setup_dir("abseq")
        d_bin = os.path.join(d, 'bin')
        plat, _ = _get_sys_info()

        if _needs_installation('clustalo'):
            b_clustal = install_clustal_omega(d)
            _syml(b_clustal, d_bin)
        else:
            print("Found clustalo, skipping installation")

        if _needs_installation('fastqc'):
            b_fastqc = install_fastqc(d)
            _syml(b_fastqc, d_bin)
        else:
            print("Found fastqc, skipping installation")

        if plat == WIN:
            if _needs_installation('flash'):
                install_flash(d)
            else:
                print("Found FLASh, skipping installation")
        else:
            if _needs_installation('leehom'):
                b_leehom = install_leehom(d)
                _syml(b_leehom, d_bin)
            else:
                print("Found leeHom, skipping installation")

        if _needs_installation('gs'):
            install_ghost_script(d)
        else:
            print("Found ghostscript, skipping installation")

        if _needs_installation('igblast'):
            retvals = install_igblast(d)
            for b in retvals:
                _syml(b, d_bin)
        else:
            print("Found igblast, skipping installation")

        # install TAMO regardless, bug fixes + custom functions / constructors used in AbSeq
        install_TAMO()

        if 'IGDATA' not in os.environ:
            with FTPBlast('ftp.ncbi.nih.gov', versions['igblast'][-1]) as blast:
                blast.download_edit_imgt_pl(d)
                igdata_dir = os.path.join(d, 'igdata')
                if not os.path.exists(igdata_dir):
                    os.makedirs(igdata_dir)
                blast.download_internal_data(igdata_dir)
                blast.download_optional_file(igdata_dir)
        else:
            print("Found IGDATA in ENV, skipping download")

        if 'IGBLASTDB' not in os.environ:
            # download human and mouse IMGT GeneDB
            download_imgt(d, "Homo+sapiens", "human")
            download_imgt(d, "Mus", "mouse")

            # create IGBLASTDB's directory
            database_dir = os.path.join(d, 'databases')
            if not os.path.exists(database_dir):
                os.makedirs(database_dir)

            # if we don't have edit_imgt_file.pl script, download it!
            if not os.path.exists(os.path.join(d, 'edit_imgt_file.pl')):
                with FTPBlast('ftp.ncbi.nih.gov', versions['igblast'][-1]) as blast:
                    blast.download_edit_imgt_pl(d)

            # if we don't have makeblastdb, download it!
            if not os.path.exists(os.path.join(d_bin, 'makeblastdb')):
                retvals = install_igblast(d)
                for b in retvals:
                    _syml(b, d_bin)

            igblast_compat(os.path.join(d, 'edit_imgt_file.pl'), os.path.join(d_bin, 'makeblastdb'),
                           os.path.join(d, 'imgt_human'),
                           os.path.join(d, 'databases'))
            igblast_compat(os.path.join(d, 'edit_imgt_file.pl'), os.path.join(d_bin, 'makeblastdb'),
                           os.path.join(d, 'imgt_mouse'), os.path.join(d, 'databases'))
        else:
            print("Found IGBLASTDB in ENV, skipping download")

        # replace install.run(self)
        # install.run(self)
        self.do_egg_install()


def readme():
    with open("README.md") as f:
        return f.read()


setup(name="AbSeq",
      version="1.1.15",
      description="Quality control pipeline for antibody libraries",
      license="placeholder",
      long_description=readme(),
      author="CSL",
      author_email="placeholder",
      maintainer="CSL",
      maintainer_email="placeholder",
      # pandas requires numpy installed, it's a known bug in setuptools - put in both setup and install requires
      # UPDATE Wed Feb 21 13:15:43 AEDT 2018 - moved into pre-installation stage
      setup_requires=['numpy>=1.11.3', 'pytz', 'python-dateutil', 'psutil', 'six'],
      install_requires=['numpy>=1.11.3', 'pandas>=0.20.1', 'biopython>=1.66', 'weblogo>=3.4', 'matplotlib>=1.5.1',
                        'tables>=3.2.3.1', 'psutil', 'matplotlib-venn', 'pyyaml'],
      packages=find_packages(),
      # NOTE TO PROGRAMMER: IF YOU CHANGE 3rd_party TO SOME OTHER DIRECTORY NAME, MAKE SURE YOU CHANGE
      # IT IN config.py AND MANIFEST.in TOO! (just search for this comment and you'll find the exact location)
      package_data={
          'abseq': ['rscripts', '3rd_party']
      },
      include_package_data=True,
      cmdclass={
          'install': ExternalDependencyInstaller,
      },
      entry_points={
          'console_scripts': ['abseq=abseq.abseqQC:main'],
      })
