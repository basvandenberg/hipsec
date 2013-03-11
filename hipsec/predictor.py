#!/usr/bin/env python

'''
Hipsec predicts if gene overexpression of an extracellular protein will lead 
to successful high-level production in Aspergillus niger. The prediction is 
based on a protein's amino acid sequence. More information about the method
can be found in:

Exploring sequence characteristics related to high-level production of secreted 
proteins in Aspergillus niger. B.A. van den Berg, M.J.T. Reinders, M. Hulsman, 
L. Wu, H.J. Pel, J.A. Roubos, D. de Ridder (2012) PLoS ONE (accepted).

Copyright (C) 2013  B.A. van den Berg

Hipsec is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>. 
'''

import argparse

__author__ = "Bastiaan van den Berg"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Bastiaan van den Berg"
__email__ = "b.a.vandenberg@tudelft.nl"

# 4 nucleotides
nucleotide_alph = 'TCAG'

# codon alphabet
codon_alph = [a+b+c for a in nucleotide_alph 
                    for b in nucleotide_alph 
                    for c in nucleotide_alph]
# remove the stop codons
codon_alph.remove('TAG')
codon_alph.remove('TAA')
codon_alph.remove('TGA')

# the 20 letter amino acid alphabet, all other letters will be ignored
aa_alph = 'ACDEFGHIKLMNPQRSTVWY'

# the amino acid weights as obtained from the trained hom classifier. These 
# values correspond to the (normalized) x-values in Figure 4A of the paper.
aa_weights = {
    'A':   1.64749, 'C': -11.99002, 'E':   1.22201, 'D':   8.39834, 
    'G':   7.33165, 'F':  16.08514, 'I':  -1.47606, 'H': -11.25671, 
    'K': -25.38251, 'M': -17.76061, 'L':  -8.52634, 'N':  20.41218, 
    'Q':   6.14166, 'P':  -4.25908, 'S':  -5.27571, 'R': -13.18086, 
    'T':   3.38485, 'W':  11.89816, 'V':   4.50645, 'Y':  18.07998
}
# the codon weights as obtained from classifier trained on ORF sequences.
codon_weights = {
    'AAA': -13.24539,  'AAC': 19.9749,   'AAG': -20.62491, 'AAT': 20.11781,
    'ACA': -12.52378,  'ACC': 9.920006,  'ACG': 2.576147,  'ACT': 1.63959,
    'AGA': -2.56792,   'AGC': -13.83788, 'AGG': -4.947894, 'AGT': -4.001029,
    'ATA': -3.565024,  'ATC': 5.284242,  'ATG': -20.41774, 'ATT': -3.060893,
    'CAA': 1.750848,   'CAC': -8.81768,  'CAG': 9.300577,  'CAT': -5.821019,
    'CCA': -8.063339,  'CCC': -6.57183,  'CCG': 5.319512,  'CCT': 2.68469,
    'CGA': -5.380213,  'CGC': -9.388379, 'CGG': -4.85336,  'CGT': -6.490896,
    'CTA': -7.571246,  'CTC': -13.8171,  'CTG': -3.836425, 'CTT': -6.602978,
    'GAA': -1.038149,  'GAC': 13.03474,  'GAG': 1.487759,  'GAT': 5.004491,
    'GCA': -0.7828718, 'GCC': -2.744247, 'GCG': 9.285551,  'GCT': -2.907307,
    'GGA': 5.406919,   'GGC': 9.610919,  'GGG': 9.563319,  'GGT': 6.019527,
    'GTA': 0.1124082,  'GTC': 6.417373,  'GTG': -1.972404, 'GTT': -4.144703,
    'TAC': 20.65884,   'TAT': 9.108252,  'TCA': -6.349592, 'TCC': 8.274822,
    'TCG': 13.55918,   'TCT': -6.69318,  'TGC': -12.96347, 'TGG': 16.10308,
    'TGT': 2.773571,   'TTA': -10.60087, 'TTC': 5.964304,  'TTG': 0.1873962,
    'TTT': 15.06293
}

# the constant as obtained from the trained hom classifier
beta_aa = -1.11158
beta_codon = -1.62544

def aa_composition(seq):
    '''
    This function returns the amino acid composition of seq. A dict with the 
    fractions for each letter in aa_alph. Only the fractions of the letters in 
    aa_alph are returned, all other letters in seq are disregarded. Therefore, 
    only if seq contains only letters from aa_alph, then the returned list of 
    fractions add up to 1.0.
    
    NOTE: the function does not check if seq is a 'correct' protein sequence!

    seq: this should be an upper case amino acid mature protein sequence 
    '''
    if(len(seq) > 0):
        comp = {}
        for letter in aa_alph:
            comp[letter] = float(seq.count(letter)) / len(seq)
        return comp
    else:
        return dict(zip(aa_alph, 20 * [0.0]))

def codon_composition(seq):
    '''
    This function returns the codon composition of seq. A dict with the 
    fractions for each codon in codon_alph. Only the fractions of the codons in 
    codon_alph are returned, all other letter triplets in seq are disregarded. 
    
    NOTE: the function does not check if seq is a 'correct' ORF sequence!

    seq: this should be an upper case nucleotide ORF sequence 
    '''
    if(len(seq) > 2):
        comp = {}
        codon_seq = window_seq(seq, 3, overlapping=False)
        for codon in codon_alph:
            comp[codon] = float(codon_seq.count(codon)) / len(codon_seq)
        return comp
    else:
        return dict(zip(codon_alph, 20 * [0.0]))

def window_seq(seq, window_size, overlapping=False):
    '''
    If the length of the sequence is not a multiple of the window size, the
    last letters are NOT returned.
    >>> s = 'AAAACCACCAAAA'
    >>> window_seq(s, 3)
    ['AAA', 'ACC', 'ACC', 'AAA']
    '''
    if(window_size < 2):
        return list(seq)
    else:
        start = 0
        stop = len(seq) - window_size + 1
        step = window_size
        if(overlapping):
            step = 1
        return [seq[i:i + window_size] for i in range(start, stop, step)]

def classification(seq, codon_sequence=False): 
    ''' 
    This method performs the classification, which is a weighted sum of 
    the sequence composition of seq plus a beta constant.
    '''
    
    if(codon_sequence):
        seq_alph = codon_alph
        weights = codon_weights
        beta = beta_codon
        comp = codon_composition(seq)
    else:
        seq_alph = aa_alph
        weights = aa_weights
        beta = beta_aa
        comp = aa_composition(seq)
    
    # calculate sum of the weighed composition
    cl = 0.0
    for letter in seq_alph:
        cl += weights[letter] * comp[letter]
    
    # add the beta constant
    return cl + beta
    
# This is a generator function
def read_fasta(handle):
    '''
    A very basic fasta parser. Not sure if it complies to the fasta standard.
    '''

    # initialize sequence id and string to an empty string
    seq_id = ''
    seq_ann = ''
    seq_str = ''

    # iterate over each line in the fasta file
    for line in handle:

        if(seq_id == '' and seq_str == ''):
            if(line[0] == '>'):
                seq_id = line.split()[0][1:]
                seq_ann = ' '.join(line.split()[1:])
            elif(line[0] == '#'):
                pass # comment
            elif(line.strip()):
                # non-empty line...
                raise(Exception, 'Error in fasta file')
        else:
            if(line == '' or line[0] == '>'):
                yield (seq_id, seq_str, seq_ann)
                seq_str = ''
                if(line[0] == '>'):
                    seq_id = line.split()[0][1:]
                else:
                    seq_id = ''
            else:
                seq_str += line.strip()

    # return the last sequence (not if the file was empty)
    if not(seq_id == ''):
        yield (seq_id, seq_str, seq_ann)

def get_output(fasta_file_n, codon_sequence=False):
    '''
    Returns a list of tuples containing the predictions of the sequences in
    the provided fasta file. The returned tuples contain three items: 
    sequence id, prediction outcome, sequence annotation.
    '''
    return [item for item in get_outputx(fasta_file_n, codon_sequence)]

# this is a generator function
def get_outputx(fasta_file_n, codon_sequence=False):
    '''
    This functions yields one prediction at the time, which could be usefull 
    in case of large fasta files. Predictions are made for the sequences in 
    the provided fasta file (fasta_file_n). The function returns tuples that
    contain three items: sequence id, prediction outcome, and sequence 
    annotation. The sequence id and sequence annotation are fetched from the 
    fasta file. The prediction outcome is calculated using the corresponding
    protein sequence.
    '''
    with open(fasta_file_n, 'r') as fin:
        for (seq_id, sequence, seq_ann) in read_fasta(fin):
            yield (seq_id, classification(sequence, codon_sequence), seq_ann)

def write_output(fasta_file_n, out_file_n, codon_sequence=False):
    '''
    This function runs the predictor on the sequences in the fasta file and 
    writes the results to the output file.
    '''
    with open(out_file_n, 'w') as fout:
        for (seq_id, cl, seq_ann) in get_outputx(fasta_file_n, codon_sequence):
            fout.write('%s\t%f\t%s\n' % (seq_id, cl, seq_ann))

def get_output_string(fasta_file_n, codon_sequence=False):
    '''
    This function runs the predictor on the sequences in the fasta file and
    returns the results as a string.
    '''
    st = ''
    for (seq_id, cl, seq_ann) in get_outputx(fasta_file_n, codon_sequence):
        st += '%s\t%f\t%s\n' % (seq_id, cl, seq_ann)
    return st

if __name__ == '__main__':
    
    # parse parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    # by default a protein amino acid sequence is expected
    parser.add_argument('-c', '--codon_sequence', action='store_true',
            default=False)
    args = parser.parse_args() 
 
    # run
    write_output(args.input_file, args.output_file, args.codon_sequence)

