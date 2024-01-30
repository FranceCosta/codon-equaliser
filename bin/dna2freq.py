#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    Converts CDS DNA to codon frequencies for a given organism. 
    Francesco Costa 28/01/2024

"""

import os
import argparse
import json
import constants
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
    "--cds",
    required=True,
    help=".fna file containing CDS of the organism of interest",
    type=pathExists,
)

parser.add_argument(
    "--output_dir",
    required=False,
    help="",
    default="data/frequency_tables",
    type=createDir,
)

parser.add_argument(
    "--organism",
    required=True,
    help="Organism to name output .json",
    type=str
)


def main():
    args = parser.parse_args()
    freq_dic = getFreqs(args.cds)
    with open(os.path.join(args.output_dir, f"{args.organism}.json"), "w") as fh:
        json.dump(freq_dic, fh)

def getFreqs(cds: str, codon_dic=constants.codons) -> dict:
    """
       
        Get codon usage frequencies. Sequences that do not start with ATG and that do not contain multiples of 3 amino 
        acids are excluded

        PARAMETERS
        ----------
        cds: str, path to file containing CDS sequences 
        codon_dic: dict, dictionary containing codons to amino acid in single letter code

        RETURNS
        -------
        dict: contains codon frequencies for each amino acid
    
    """
    aa_dic = {aa:0 for aa in set(codon_dic.values())}
    freq_dic = {codon: 0  for codon in codon_dic}
    sequences = list(SeqIO.parse(cds, "fasta"))
    for sequence in sequences:
        seq = str(sequence.seq)
        # check that the sequence starts with start codon and has multiples of 3
        cond1 = seq[:3] == "ATG"
        cond2 = len(seq) % 3 == 0
        if cond1 and cond2:
            for index in range(0, len(seq), 3):
                codon = seq[index:index+3]
                # exclude codons with non canonical characters (N is seldomly present)
                if codon in codon_dic.keys():
                    freq_dic[codon] += 1
    
    # now convert this into frequencies
    for codon in freq_dic:
        aa = codon_dic[codon]
        aa_dic[aa] += freq_dic[codon]
    
    new_freq_dic = {aa: {} for aa in set(codon_dic.values())}
    for codon in freq_dic:
        aa = codon_dic[codon] # find the amino acid
        total = aa_dic[aa] # find total
        freq = freq_dic[codon] / total
        new_freq_dic[aa][codon] = round(freq, 3)
    
    return new_freq_dic



if __name__ == "__main__":
    main()