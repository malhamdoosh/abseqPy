import sys
import os

def fastq(it):
    seq_id = None
    try:
        while 1:
            seq_id = it.next().rstrip()
            assert seq_id[0] == '@'

            seq_buf = []

            while 1:
                seq = it.next().rstrip()
                if seq.startswith('+'):
                    break
                seq_buf.append(seq)

            if len(seq) != 1:
                assert seq[1:] == seq_id[1:]

            seq = ''.join(seq_buf)

            qual_buf = []
            qual_count = 0
            while qual_count < len(seq):
                qual = it.next().rstrip()
                qual_buf.append(qual)
                qual_count += len(qual)

            assert qual_count == len(seq)

            qual = ''.join(qual_buf)

            yield seq_id, seq, qual
            seq_id = None
    except StopIteration:
        if seq_id:
            raise exceptions.RuntimeError('partial fastq sequence')



def fasta(it):
    seq_id = None
    seq_buf = []

    try:
        seq_id = it.next().rstrip()
        while 1:
            assert seq_id[0] == '>'

            seq_buf = []

            while 1:
                seq = it.next().rstrip()
                if seq.startswith('>'):
                    next_seq_id = seq
                    break
                seq_buf.append(seq)

            yield (seq_id, ''.join(seq_buf))

            seq_id = next_seq_id
    except StopIteration:
        if len(seq_buf):
            yield (seq_id, ''.join(seq_buf))



__all__ = [
    'fastq',
    'fasta'
]
