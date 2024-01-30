#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    Equalise sequence using closest codons in frequency to host organism
    Francesco Costa 29/01/2024

"""


import os
import argparse
import json
from constants import codons
from Bio import SeqIO

def pathExists(path: str) -> str:
    """Validate input path"""
    if os.path.exists(path):
        return path
    else:
        raise NotADirectoryError(path)


def createDir(dirPath: str) -> str:
    """Cretes directory if non existing"""
    os.makedirs(dirPath, exist_ok=True)
    return dirPath


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "--sequences",
    required=True,
    help=".fasta file containing sequences to equalise",
    type=pathExists,
)

parser.add_argument(
    "--origin_organism",
    required=False,
    help="Organism of origin",
    default="hsapiens",
    type=str,
)

parser.add_argument(
    "--dest_organism",
    required=False,
    default="ecoli",
    help="Organism of destination",
    type=str,
)

parser.add_argument(
    "--freq_tables",
    required=False,
    help=".json tables with codon frequencies",
    default="data/frequency_tables",
    type=pathExists,
)


def main():
    args = parser.parse_args()
    with open(os.path.join(args.freq_tables, f"{args.origin_organism}.json"), 'r') as fh:
       origin_dic = json.load(fh)
    with open(os.path.join(args.freq_tables, f"{args.dest_organism}.json"), 'r') as fh:
        dest_dic = json.load(fh)
    sequences = list(SeqIO.parse(args.sequences, "fasta"))
    
    output_sequences = []
    for sequence in sequences:
        out_sequence = equaliser(sequence, origin_dic, dest_dic)
        output_sequences.append(out_sequence)
    # run equaliser

def equaliser(sequence:SeqIO.SeqRecord, origin_dic:dict, dest_dic:dict, codon_table=codons) -> SeqIO.SeqRecord:
    """
    
        Return a sequence with closer composition to the original as possible.

        PARAMETES
        --------
        sequence: 
        origin_dic:
        dest_dic:
        codon_table:

        RETURNS
        -------
        sequence:
    
    """
    assert len(sequence.seq) % 3 == 0, "Sequence lenght should be multiple of 3"
    assert any([sequence.seq[i:i+3] not in codon_table.keys() for i in range(0, len(sequence.seq), 3)]) == False,\
    "Non-canonical codon detected"
    
    origin_freqs = []
    dest_freqs = []
    output_sequence = []
    
    for indx in range(0, len(sequence.seq), 3):
        codon = sequence.seq[indx:indx+3]
        aa = codon_table[codon]
        origin_freq = origin_dic[aa][codon]
        # find the closest destination freq
        diff = 1
        for codon in dest_dic[aa]:
            new_diff = abs(origin_freq - dest_dic[aa][codon])
            if new_diff < diff:
                diff = new_diff
                dest_codon = codon
                dest_freq = dest_dic[aa][codon]
        print(aa, origin_freq, dest_freq)

        

if __name__ == "__main__":
    main()