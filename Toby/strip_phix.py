#!/usr/bin/env python
import sys
import os
import io
import argparse
import string
import itertools

from seqfile import *
from argparse_utils import *

def revcomp(seq):
    return string.translate(seq[::-1], string.maketrans('XACMGRSVTWYHKDBNxacmgrsvtwyhkdbn', 'XTGKCYSBAWRDMHVNxtgkcysbawrdmhvn'))

def make_kmer_set(seq, k):
    out = set()
    for i in xrange(len(seq)-k):
        out.add(seq[i:i+k])
    return out

if __name__ == '__main__':
    parser = argparse.ArgumentParser('filter PhiX sequences from input fastq files')
    parser.add_argument('--phix',   '-p', type = argparse.FileType('r'),   help='PhiX fasta file',        default=open('phix/phix.fa'))
    parser.add_argument('--kmer',   '-k', type = int,                      help='kmer size for matching', default=8)
    parser.add_argument('--reads1', '-1', type = MaybeCompressedFile('r'), help='first fastq read pair',  nargs='*')
    parser.add_argument('--reads2', '-2', type = MaybeCompressedFile('r'), help='second fastq read pair', nargs='*')
    parser.add_argument('--reads',  '-r', type = MaybeCompressedFile('r'), help='unpaired fastq',         nargs='*')

    args = parser.parse_args()

    if len(args.reads1) != len(args.reads2):
        print >>sys.stderr, 'paired fastq inputs are unbalanced'
        sys.exit(1)

    phix_id, phix_seq = fasta(args.phix).next()

    phix_kmers = make_kmer_set(phix_seq, args.kmer) | make_kmer_set(revcomp(phix_seq), args.kmer)

    for r1, r2 in zip(args.reads1, args.reads2):
        for seq1, seq2 in itertools.izip(fastq(r1), fastq(r2)):
            k = make_kmer_set(seq1[1], args.kmer) | make_kmer_set(seq2[1], args.kmer)
            print seq1[0], len(k & phix_kmers)
