"""
PermuteTools.py -- Utility functions for enumerationg lists of sequence permutations

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
"""

from TAMO.MotifTools import revcomplement as RC

def permute(letters, depth, seqs=[''],curdepth=0):
    """
    permute(letters, depth, seqs=[''],curdepth=0) -- Generate all possible sequences from the alphabet
                                                     "letters" of length "depth" (sorry about the
                                                     variable names).  For example, "permute('ACGT',4)"
                                                     generates all 256 possible 4-letter DNA sequences.
    """
    newseqs = []
    for seq in seqs:
        for letter in letters:
            newseqs.append(seq + letter)
    if depth > curdepth:
        return(permute(letters,depth,newseqs,curdepth + 1))
    else:
        return(seqs)

	
def restricted_permute(base_letters, extra_letters, depth, seqs = ['']):
    """
    restricted_permute(base_letters, extra_letters, depth, seqs = ['']) -- [Utility function] Like permute()
                                                                        above, but allows at most one letter
                                                                        from "extra_letters" into the
                                                                        permutations.  Used for including
                                                                        ambiguity codes. 
    """
    if depth == 0: return seqs
    ans = []
    extended = base_letters + extra_letters
    for seq in seqs:
        used_up = 0
        for L in extra_letters:
            if L in seq: used_up = 1
        if used_up: Ls = base_letters
        else:       Ls = extended
        for L in Ls:
            ans.append(seq + L)
    return restricted_permute(base_letters, extra_letters, depth-1, ans)
        
def uniq_syl_pairs(width):
    """
    uniq_syl_pairs(width) -- Generates all possible (left,right) sequences, allowing at most one
                             ambiguous letter (S, W, R, or Y) in either left or right.  uniq_syl_pairs(4)
                             returns 836,160 tuples, beginning ('AAAA', 'AAAA'), ('AAAA', 'AAAC'),
                             ('AAAA', 'AAAG'), ('AAAA', 'AAAT'), ('AAAA', 'AAAS'), ...
    """
    syls = restricted_permute(list('ACGT'), list('SWRY'),width)
    memo = {}
    tups = []
    for left in syls:
        for right in syls:
            word = left+right
            rc   = RC(word)
            if memo.has_key(word) or memo.has_key(rc): continue
            memo[word] = 1
            tups.append((left,right))
    return tups

def dimer_words(width):
    """
    dimer_words(width) -- For every possible word of width "width" with at most one ambiguity code,
                          assemble the words of form word-gap-(word'), where (word') is the reverse
                          complement of word, and gap is a string of "N"s ranging in length from
                          0 to 12.
                          
    """
    syls = restricted_permute(list('ACGT'), list('SWRY'),width)
    words = []
    memo  = {}
    tups  = []
    for syl in syls:
        for af in [syl, RC(syl)]:
            word = syl+af
            if memo.has_key(word) or memo.has_key(RC(word)): continue
            memo[word] = 1
            tups.append((syl,af))
    for gaplen in range(12):
        gap = 'N'*gaplen
        for left, right in tups:
            words.append('%s%s%s'%(left,gap,right))
    return words

def gapped_words(width):
    """
    gapped_words(width) -- For every pair returned from uniq_syl_pairs (see documentation for that function)
                           assemble into words with spacings 0-12.   So (AA, GT) --> AAGT, AANGT, AANNGT, ...
    """
    tups = uniq_syl_pairs(width)
    words = []
    for gaplen in range(12):
        gap = 'N'*gaplen
        print gap, len(words)
        for left, right in tups:
            words.append('%s%s%s'%(left,gap,right))
    return words

def dimer_groups(width):
    """
    dimer_groups(width) -- For every possible word of width "width" with at most one ambiguity code,
                           assemble the REGULAR EXPRESSIONS of form word-gap-(word'), where (word')
                           is the reverse complement of word, and the gap is the text string '.{a,b}'
                           where 'a' and 'b' vary between 0-16 such that b > a.  Intended to represent
                           all possible variable gap motifs that are symmetric.
    """
    syls = restricted_permute(list('ACGT'), list('SWRY'),width)
    words = []
    memo  = {}
    tups  = []
    for syl in syls:
        for af in [syl, RC(syl)]:
            word = syl+af
            if memo.has_key(word) or memo.has_key(RC(word)): continue
            memo[word] = 1
            tups.append((syl,af))
    for startgap in range(15):
        for stopgap in range(startgap,16):
            for left, right in tups:
                words.append('%s.{%d,%d}%s'%(left,startgap,stopgap,right))
    return words

def expt():
    words = gapped_words(3)
    print len(words)
    
        

            
                
        
    
