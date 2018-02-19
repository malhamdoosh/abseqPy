from __future__ import print_function

import sys
import glob
import os
import platform
import zipfile
import struct
import re
from Bio import SeqIO

from subprocess import check_output
from ftplib import FTP

from config import EXTERNAL_DEP_DIR

MAC = 'Darwin'
LIN = 'Linux'
WIN = 'Windows'

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
        installation_path = (installation_dir + '/' + binary).replace('//', '/')
        with open(installation_path, "wb") as fp:
            self.ftp.retrbinary('RETR ' + path, fp.write)

        old_dir = os.path.abspath(".")
        os.chdir(installation_dir)
        _ = check_output(['tar', '-xvzf', binary])
        os.chdir(old_dir)

        return glob.glob(installation_dir + '/ncbi-igblast-' + self.version + '/bin/*')

    def download_edit_imgt_pl(self, download_dir):
        path = '/blast/executables/igblast/release/edit_imgt_file.pl'
        download_path = (download_dir + '/edit_imgt_file.pl').replace('//', '/')
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
            download_path = (download_dir + '/internal_data/' + s + '/').replace('//', '/')
            os.makedirs(download_path)
            for filename in filenames:
                # ignore rhesus_monkey's CVS directory
                if filename == 'CVS':
                    continue
                with open(download_path + filename, "wb") as fp:
                    self.ftp.retrbinary('RETR ' + filename, fp.write)
            self.ftp.cwd('../')


    def download_optional_file(self, download_dir):
        path = '/blast/executables/igblast/release/optional_file/'
        self.ftp.cwd(path)
        filenames = self.ftp.nlst()
        download_path = (download_dir + '/optional_file/').replace('//', '/')
        os.makedirs(download_path)
        for filename in filenames:
            with open(download_path + filename, "wb") as fp:
                self.ftp.retrbinary('RETR ' + filename, fp.write)


def _get_sys_info():
    return platform.system(), 8 * struct.calcsize("P")

def _error(msg, stream=sys.stderr, abort=1):
    print(msg, file=stream)
    if abort:
        sys.exit(abort)

def _syml(src, dest):
    binary_name = os.path.basename(src)
    if src:
        os.symlink(os.path.abspath(src), (dest + '/' + binary_name).replace('//', '/'))

def setup_dir(root):
    output = (root + "/" + EXTERNAL_DEP_DIR).replace('//', '/')
    if not os.path.exists(output):
        os.makedirs(output)
    return output

def install_clustal_omega(installation_dir="."):
    plat, bit = _get_sys_info()

    # clustalo needs to create a dir
    installation_dir = (installation_dir + '/' + 'clustal-omega').replace('//', '/')
    if not os.path.exists(installation_dir):
        os.makedirs(installation_dir)

    if plat == MAC:
        addr = 'http://www.clustal.org/omega/clustal-omega-1.2.3-macosx'
    elif plat == LIN:
        if bit == 64:
            addr = 'http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64'
        elif bit == 32:
            addr = 'http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-32-bit'
        else:
            _error('Unknown architecture. Detected a non 32 or 64 bit system.')
    elif plat == WIN:
        addr = 'http://www.clustal.org/omega/clustal-omega-1.2.2-win64.zip'
    else:
        _error('Unknown system architecture. Non windows, mac or linux detected')
    binary = (installation_dir + '/' + 'clustalo').replace('//', '/')
    # install binary
    _ = check_output(['curl', addr, '-o', binary])
    # add execution bit
    os.chmod(binary, 0o777)
    return binary


def install_fastqc(installation_dir="."):
    addr = 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip'
    zipname = (installation_dir + '/' + addr.split('/')[-1].strip()).replace('//', '/')
    _ = check_output(['curl', addr, '-o', zipname])
    unzipped_name = 'FastQC'
    zip_ref = zipfile.ZipFile(zipname, 'r')
    zip_ref.extractall(installation_dir)
    zip_ref.close()
    binary = (installation_dir + '/' + unzipped_name + '/' + 'fastqc').replace('//', '/')
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

    return (installation_dir + '/leeHom' + '/src/leeHomMulti').replace('//', '/')
    

def install_ghost_script(installation_dir='.', threads=2):
    addr = 'https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs922/ghostpdl-9.22.tar.gz'
    tarname = addr.split('/')[-1].strip()
    old_dir = os.path.abspath('.')

    target_dir = os.path.abspath(installation_dir)

    os.chdir(installation_dir)
    _ = check_output(['curl', '-L', addr, '-o', tarname])
    _ = check_output(['tar', '-xvzf', tarname])
    ghs_dir = os.path.splitext(os.path.splitext(tarname)[0])[0]
    os.chdir(ghs_dir)
    _ = check_output(['./configure', '--prefix={}'.format(target_dir)]) 
    _ = check_output(['make', '-j', str(threads)])
    _ = check_output(['make', 'install'])
    os.chdir(old_dir)

    # dont need to return binary directory, it's already in installation_dir/bin


def install_igblast(installation_dir='.', version='1.8.0'):
    plat, bit = _get_sys_info()

    with FTPBlast('ftp.ncbi.nih.gov', version) as blast:
        if plat == MAC:
            bins = blast.install_bins('ncbi-igblast-{}-x64-macosx.tar.gz'.format(version), installation_dir)
        elif plat == WIN:
            bins = blast.install_bins('ncbi-igblast-{}-x64-win64.tar.gz'.format(version), installation_dir)
        elif plat == LIN:
            bins = blast.install_bins('ncbi-igblast-{}-x64-linux.tar.gz'.format(version), installation_dir)
        else:
            _error("Unknown platform detected")

        blast.download_edit_imgt_pl(installation_dir)
        igdata_dir = (installation_dir + '/igdata').replace('//', '/')
        os.makedirs(igdata_dir)
        blast.download_internal_data(igdata_dir)
        blast.download_optional_file(igdata_dir)

    return bins


def download_imgt(download_dir, species, species_layman):
    links = [
            "http://www.imgt.org/genedb/GENElect?query=7.14+IGHV&species=",
            "http://www.imgt.org/genedb/GENElect?query=7.14+IGHD&species=",
            "http://www.imgt.org/genedb/GENElect?query=7.14+IGHJ&species=",
            "http://www.imgt.org/genedb/GENElect?query=7.14+IGKV&species=",
            "http://www.imgt.org/genedb/GENElect?query=7.14+IGKJ&species=",
            "http://www.imgt.org/genedb/GENElect?query=7.14+IGLV&species=",
            "http://www.imgt.org/genedb/GENElect?query=7.14+IGLJ&species="
    ]

    path = (download_dir + '/imgt_' + species_layman + "/").replace('//', '/')
    os.makedirs(path)
    for url in links:
        gene = url[url.find("+") + 1:url.find("&")].lower()
        output = "{}_{}.imgt.raw".format(path + species_layman, gene)
        _ = check_output(['curl', '-L', url + species, '-o', output])
        # TODO: parse file to get pure genes only
        with open(output[:output.rfind(".")], "w") as writer, \
             open(output) as reader:

            line = reader.readline()
            while not line.startswith("<b>Number of results"):
                line = reader.readline()
            if not line:
                raise Exception("File has no IMGT sequences")

            reader.readline() #\n
            reader.readline() #<pre>

            for line in reader:
                if line.startswith("</pre>"):
                    # finish writing sequences
                    break
                writer.write(line)
            # remove raw
            os.remove(output)
            
def igblast_compat(edit_imgt_bin, make_blast_bin, data_dir, output_dir):
    for f in os.listdir(data_dir):
        clean_fasta = output_dir+ 'imgt_' + f[:f.find(".")]
        os.system(edit_imgt_bin+ ' ' + data_dir + f + ' > ' + clean_fasta)
        records = []
        for rec in SeqIO.parse(clean_fasta, 'fasta'):
            rec.description = ''
            rec.seq = rec.seq.upper()
            records.append(rec)
        SeqIO.write(records, clean_fasta, 'fasta')
        _ = check_output([make_blast_bin, '-parse_seqids', '-dbtype', 'nucl', '-in', clean_fasta])
        if len(re.findall('ig[h,k,l][v,c]', clean_fasta)) > 0:
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
        


if __name__ == '__main__':

    # create external deps dir
    d = setup_dir(".")
    d_bin = (d + '/bin').replace('//', '/')
    os.makedirs(d_bin)
    
    # BIN: d + 'clustal-omega/clustalo'
    b_clustal = install_clustal_omega(d)
    _syml(b_clustal, d_bin)

    # BIN: d + 'FastQC/fastqc'
    b_fastqc = install_fastqc(d)
    _syml(b_fastqc, d_bin)

    # BIN: d + 'leeHom/src/leeHomMulti'
    b_leehom = install_leehom(d)
    _syml(b_leehom, d_bin)

    install_ghost_script(d)
    retvals = install_igblast(d)
    for b in retvals:
        _syml(b, d_bin)

    # download human and mouse IMGT GeneDB
    download_imgt(d, "Homo+sapiens", "human")
    download_imgt(d, "Mus", "mouse")

    # create IGBLASTDB's directory
    os.makedirs(d + '/databases/')
    igblast_compat(d + '/edit_imgt_file.pl', d_bin + '/makeblastdb', d + '/imgt_human/', d + '/databases/')
    igblast_compat(d + '/edit_imgt_file.pl', d_bin + '/makeblastdb', d + '/imgt_mouse/', d + '/databases/')
