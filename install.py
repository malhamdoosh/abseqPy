from __future__ import print_function

import os
import sys
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


def _save_as(url, fname, chmod=True):
    import requests
    r = requests.get(url)
    attempts = 0

    # keep trying until we get it
    while r.status_code != 200 and attempts < 10:
        r = requests.get(url)
        attempts += 1

    if r.status_code != 200:
        raise Exception("Cannot download {}, fatal error".format(url))

    with open(fname, 'wb') as fp:
        for c in r:
            fp.write(c)
    if chmod:
        os.chmod(fname, 0o777)


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
    if isinstance(v, bool) or isinstance(software_version, bool):
        return software_version != v
    if isinstance(v, list):
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


def _setup_dir(root):
    root = os.path.abspath(root)
    if not os.path.exists(root):
        os.makedirs(root)
    return root


def _install_clustal_omega(installation_dir=".", version=versions['clustalo'][-1]):
    # can't use versions yet, pre-compiled binaries are a little out of sync

    plat, bit = _get_sys_info()
    # clustalo needs to create a dir
    clustalo_installation_dir = os.path.join(installation_dir, 'clustal-omega')
    if not os.path.exists(clustalo_installation_dir):
        os.makedirs(clustalo_installation_dir)

    binary = os.path.join(clustalo_installation_dir, 'clustalo')
    if plat == MAC:
        addr = 'http://www.clustal.org/omega/clustal-omega-1.2.3-macosx'
        _save_as(addr, binary)
    elif plat == LIN:
        if bit == 64:
            addr = 'http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64'
        elif bit == 32:
            addr = 'http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-32-bit'
        else:
            _error('Unknown architecture. Detected a non 32 or 64 bit system.')
        # noinspection PyUnboundLocalVariable
        _save_as(addr, binary)
    elif plat == WIN:
        windows_bin = 'clustal-omega-1.2.2-win64.zip'
        addr = 'http://www.clustal.org/omega/' + windows_bin
        _save_as(addr, windows_bin, chmod=False)
        zip_ref = zipfile.ZipFile(windows_bin)
        zip_ref.extractall(clustalo_installation_dir)
        zip_ref.close()
        # windows has no symlink, go straight to bin directory!
        for f in os.listdir(os.path.join(clustalo_installation_dir, windows_bin[:windows_bin.find('.zip')])):
            src_ = os.path.join(clustalo_installation_dir, windows_bin[:windows_bin.find('.zip')], f)
            shutil.move(src_, os.path.join(installation_dir, 'bin', f))
    else:
        _error('Unknown system architecture. Non windows, mac or linux detected')

    return binary


def _install_fastqc(installation_dir=".", version=versions['fastqc'][-1]):
    plat, _ = _get_sys_info()
    addr = 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v{}.zip'.format(version)
    zipname = os.path.join(installation_dir, os.path.basename(addr).strip())
    _save_as(addr, zipname, chmod=False)
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


def _install_leehom(installation_dir='.'):
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


def _install_flash(installation_dir='.'):
    addr = "http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11-windows-bin.zip"
    flash_ins_dir = os.path.join(installation_dir, "flash")
    if not os.path.exists(flash_ins_dir):
        os.makedirs(flash_ins_dir)
    flash_zip = 'flash.zip'
    _save_as(addr, flash_zip, chmod=False)
    zip_ref = zipfile.ZipFile(flash_zip)
    zip_ref.extractall(flash_ins_dir)
    for f in os.listdir(flash_ins_dir):
        shutil.move(os.path.join(flash_ins_dir, f), os.path.join(installation_dir, 'bin', f))


def _install_ghost_script(installation_dir='.', threads=2, version=versions['gs'][-1]):
    plat, bit = _get_sys_info()
    target_dir = os.path.abspath(installation_dir)

    if plat != WIN:
        addr = 'https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs{}/ghostpdl-{}.tar.gz'.format(
            version.replace('.', ''), version)
        tarname = os.path.basename(addr)
        old_dir = os.path.abspath('.')

        os.chdir(installation_dir)
        _save_as(addr, tarname, chmod=False)
        _ = check_output(['tar', '-xvzf', tarname])
        ghs_dir = os.path.splitext(os.path.splitext(tarname)[0])[0]
        os.chdir(ghs_dir)
        _ = check_output(['./configure', '--prefix={}'.format(target_dir)])
        _ = check_output(['make', '-j', str(threads)])
        _ = check_output(['make', 'install'])
        os.chdir(old_dir)
    # 1/07/2018 - Window's ghostscript renderer has issues removed from abseqPy's diversity_analysis::weblogo
    # else:
        # addr = "http://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs{}/ghostpcl-{}-win64.zip" \
        #     .format(version.replace('.', ''), version)
        # gs_dir = "ghostpcl-{}-win64".format(version)
        # _save_as(addr, gs_dir + ".zip", chmod=False)
        # zip_ref = zipfile.ZipFile(gs_dir + ".zip")
        # zip_ref.extractall(gs_dir)
        # # nested folder, not a typo
        # gs_out_dir = os.path.join(gs_dir, gs_dir)
        # for f in os.listdir(gs_out_dir):
        #     shutil.move(os.path.join(gs_out_dir, f), os.path.join(installation_dir, 'bin', f))
    # dont need to return binary directory, it's already in installation_dir/bin


def _install_igblast(installation_dir='.', version=versions['igblast'][-1]):
    plat, _ = _get_sys_info()

    with FTPBlast('ftp.ncbi.nih.gov', version) as blast:
        if plat == MAC:
            bins = blast.install_bins('ncbi-igblast-{}-x64-macosx.tar.gz'.format(version), installation_dir)
        elif plat == WIN:
            bins = blast.install_bins('ncbi-igblast-{}-x64-win64.tar.gz'.format(version), installation_dir)
            for bin_ in bins:
                shutil.move(bin_, os.path.join(installation_dir, 'bin', os.path.basename(bin_)))
        elif plat == LIN:
            bins = blast.install_bins('ncbi-igblast-{}-x64-linux.tar.gz'.format(version), installation_dir)
        else:
            _error("Unknown platform detected")
    # noinspection PyUnboundLocalVariable
    return bins


# no longer used. TAMO is installed by default together with abseq's setup.py
def _install_TAMO():
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


def _download_imgt(download_dir, species, species_layman):
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
        _save_as(url + species, output, chmod=False)
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


def _igblast_compat(edit_imgt_bin, make_blast_bin, data_dir, output_dir):
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


def install(directory):
    igdata_downloaded = False
    igblastdb_downloaded = False

    # create external deps dir
    if len(sys.argv) != 2:
        _error("ERROR: Usage: python {} install_path".format(sys.argv[0]))

    # if not os.path.exists('TAMO.tar.gz'):
    #     _error("ERROR: Please use this script in the directory where TAMO.tar.gz is in")

    # although setup() has this, it's installed locally in abseq's installation dir.
    # By pip.installing here, it's going to be available globally
    setup_requires = ['requests', 'numpy>=1.11.3', 'biopython>=1.66']
    for pack in setup_requires:
        # pip.main(['install', pack]) no longer supported in pip >= 10
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', pack])

    d = _setup_dir(directory)

    d_bin = os.path.join(d, 'bin')
    if not os.path.exists(d_bin):
        os.makedirs(d_bin)
        
    plat, _ = _get_sys_info()

    if _needs_installation('clustalo'):
        b_clustal = _install_clustal_omega(d)
        _syml(b_clustal, d_bin)
    else:
        print("Found clustalo, skipping installation")

    fastqc_downloaded = False
    if _needs_installation('fastqc'):
        b_fastqc = _install_fastqc(d)
        _syml(b_fastqc, d_bin)
        fastqc_downloaded = True
    else:
        print("Found fastqc, skipping installation")

    if plat == WIN:
        if _needs_installation('flash'):
            _install_flash(d)
        else:
            print("Found FLASh, skipping installation")
    else:
        if _needs_installation('leehom'):
            b_leehom = _install_leehom(d)
            _syml(b_leehom, d_bin)
        else:
            print("Found leeHom, skipping installation")

    if _needs_installation('gs'):
        _install_ghost_script(d)
    else:
        print("Found ghostscript, skipping installation")

    if _needs_installation('igblast'):
        retvals = _install_igblast(d)
        for b in retvals:
            _syml(b, d_bin)
    else:
        print("Found igblast, skipping installation")

    # install TAMO regardless, bug fixes + custom functions / constructors used in AbSeq
    # _install_TAMO()

    with FTPBlast('ftp.ncbi.nih.gov', versions['igblast'][-1]) as blast:
        if "IGDATA" in os.environ:
            igdata_path_contents = os.listdir(os.path.expandvars("$IGDATA"))
            duplicate_igdata = True
        else:
            igdata_path_contents = []
            duplicate_igdata = False

        igdata_dir = os.path.join(d, 'igdata')
        if 'internal_data' not in igdata_path_contents:
            if not os.path.exists(igdata_dir):
                os.makedirs(igdata_dir)
            blast.download_internal_data(igdata_dir)
            igdata_downloaded = True

        if 'optional_file' not in igdata_path_contents:
            if not os.path.exists(igdata_dir):
                os.makedirs(igdata_dir)
            blast.download_optional_file(igdata_dir)
            igdata_downloaded = True

    if 'IGBLASTDB' not in os.environ:
        # download human and mouse IMGT GeneDB
        _download_imgt(d, "Homo+sapiens", "human")
        _download_imgt(d, "Mus", "mouse")

        # create IGBLASTDB's directory
        database_dir = os.path.join(d, 'databases')
        if not os.path.exists(database_dir):
            os.makedirs(database_dir)

        # if we don't have edit_imgt_file.pl script, download it!
        if not os.path.exists(os.path.join(d, 'edit_imgt_file.pl')):
            with FTPBlast('ftp.ncbi.nih.gov', versions['igblast'][-1]) as blast:
                blast.download_edit_imgt_pl(d)

        # if we don't have makeblastdb, download it!
        if not os.path.exists(os.path.join(d_bin, _binary_file('makeblastdb'))):
            retvals = _install_igblast(d)
            for b in retvals:
                _syml(b, d_bin)

        _igblast_compat(os.path.join(d, 'edit_imgt_file.pl'), os.path.join(d_bin, 'makeblastdb'),
                        os.path.join(d, 'imgt_human'),
                        os.path.join(d, 'databases'))
        _igblast_compat(os.path.join(d, 'edit_imgt_file.pl'), os.path.join(d_bin, 'makeblastdb'),
                        os.path.join(d, 'imgt_mouse'), os.path.join(d, 'databases'))
        igblastdb_downloaded = True
    else:
        print("Found IGBLASTDB in ENV, skipping download")

    return fastqc_downloaded, igblastdb_downloaded, igdata_downloaded, duplicate_igdata


def ask_permission(directory):
    msg = "\nYou are about to install abseqPy's external dependencies in {}.\nIt is your \
responsibility to make sure that the folder is empty or \nthat the downloaded dependencies \
will not override your existing files.\n\nProceed? [y/N]: ".format(directory)
    try:
        inputc = raw_input
    except NameError:
        inputc = input
    try:
        return str(inputc(msg)).lower() in ['y', 'yes']
    except KeyboardInterrupt:
        return False


def _binary_file(binary):
    plat, _ = _get_sys_info()
    if plat == WIN:
        return binary + ".exe"
    return binary


def main():
    if len(sys.argv) != 2:
        print("Usage: {} <installation directory>".format(sys.argv[0]), file=sys.stderr)
        return 0

    directory = sys.argv[1]
    if ' ' in directory:
        print("Installation directory is not allowed to contain spaces!")
        return 1

    proceed = ask_permission(directory)

    if proceed:

        fastqc_downloaded, igblastdb_downloaded, igdata_downloaded, duplicate_igdata = install(directory)

        print("", file=sys.stderr)
        print("Installation complete, remember to add the following line(s) to your ~/.bashrc or equivalent",
              file=sys.stderr)
        print("", file=sys.stderr)
        if igblastdb_downloaded:
            print("\texport IGBLASTDB=\"{}\"".format(os.path.join(os.path.abspath(directory), "databases")),
                  file=sys.stderr)
        if igdata_downloaded:
            if duplicate_igdata:
                print(
                    "\nDetected existing IGDATA=\"{ori}\" in your environment.\nHowever, you did not have all "
                    "the required files in that directory.\nThe required files are downloaded in {dup},"
                    "\nyou should move these file(s) to {ori} in order for abseqPy to work.\n"
                    "\nIf you do not have permission to do so, you can copy all the files from\n{ori} to {dup},"
                    " and append\n\n\texport IGDATA=\"{dup}\"\n\nto your ~/.bashrc (or equivalent) instead.".format(
                        ori=os.path.expandvars("$IGDATA"), dup=os.path.join(os.path.abspath(directory), "igdata"))
                    , file=sys.stderr)
            else:
                print("\texport IGDATA=\"{}\"".format(os.path.join(os.path.abspath(directory), "igdata")),
                      file=sys.stderr)
        if platform.system() == "Windows":
            # windows needs extra FASTQCROOT export - perl script needs to be invoked manually, not via shebang
            if fastqc_downloaded:
                print("\texport FASTQCROOT=\"{}\"".format(os.path.join(os.path.abspath(directory), "bin")))

            # export now doesn't show ':' as Unix systems does - don't wanna confuse reader
            print("\texport PATH=\"{}\"".format(os.path.join(os.path.abspath(directory), "bin")), file=sys.stderr)
            print("\nFor windows users, export <Variable>=<Value> means you should type key=value pair\n"
                  "into your user environment variable. This can be accessed via\n\n"
                  "\tStartMenu > typing in 'env' in the search box > Click on 'Edit the system environment variables' "
                  ">\n"
                  "\tClick on 'Environment Variables...' button > Click 'New...' to add new Variable(all but PATH)\n"
                  "\tor 'Edit...' to append your PATH variable (remember to use ';'). Ask your local administrator\n"
                  "\tif this message is confusing to you.")
        else:
            print("\texport PATH=\"${{PATH}}:{}\"".format(os.path.join(os.path.abspath(directory), "bin")),
                  file=sys.stderr)

        print("", file=sys.stderr)
    else:
        print("Aborted installer.py", file=sys.stderr)

    return 0


if __name__ == '__main__':
    sys.exit(main())
