#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    Equalise sequence using closest codons in frequency to host organism
    Francesco Costa 29/01/2024

"""


import os
import argparse
import json
from Bio import SeqIO

codons = {
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TTC': 'F',
    'TTT': 'F',
    'TTA': 'L',
    'TTG': 'L',
    'TAC': 'Y',
    'TAT': 'Y',
    'TAA': '*',
    'TAG': '*',
    'TGC': 'C',
    'TGT': 'C',
    'TGA': '*',
    'TGG': 'W',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CAC': 'H',
    'CAT': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'ATA': 'I',
    'ATC': 'I',
    'ATT': 'I',
    'ATG': 'M',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'AAC': 'N',
    'AAT': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'AGC': 'S',
    'AGT': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'GAC': 'D',
    'GAT': 'D',
    'GAA': 'E',
    'GAG': 'E', 
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G'
}

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
    
    output_data = []
    for sequence in sequences:
        out_ = equaliser(sequence, origin_dic, dest_dic)
        output_data.append(out_)
    
    output_id = args.sequences.split(".")[0]
    output_sequences = [i[0] for i in output_data]
    with open(f"{output_id}_equalised.fa", "w") as fh:
        for seq in output_sequences:
            fh.write(f">{seq.id}\n{seq.seq}\n")
    
    

def equaliser(sequence:SeqIO.SeqRecord, origin_dic:dict, dest_dic:dict, codon_table=codons) -> list:
    """
    
        Return a sequence with closer composition to the original as possible.

        PARAMETES
        --------
        sequence: 
        origin_dic: codon usage dictionary of origin organism
        dest_dic: codon usage dictionary of destination organism
        codon_table: table of codons

        RETURNS
        -------
        [sequence: SeqIO.SeqRecord, origin_freqs: list, dest_freqs: list]
    
    """
    assert len(sequence.seq) % 3 == 0, "Sequence lenght should be multiple of 3"
    assert any([sequence.seq[i:i+3] not in codon_table.keys() for i in range(0, len(sequence.seq), 3)]) == False,\
    "Non-canonical codon detected"
    
    origin_freqs = []
    dest_freqs = []
    output_sequence = ""
    
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
        origin_freqs.append(origin_freq)
        dest_freqs.append(dest_freq)
        output_sequence += dest_codon
        #print(aa, origin_freq, codon, dest_freq, dest_codon)
    
    return [SeqIO.SeqRecord(output_sequence, sequence.id+"_equalised"), origin_freqs, dest_freqs]

        

if __name__ == "__main__":
    main()