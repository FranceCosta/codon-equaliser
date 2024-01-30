# Web tool to optimise protein sequences for heterologous expression

Protein expression success strongly depends on translation speed [link](https://doi.org/10.1016/j.jmb.2023.168384). 
So far, heterologous protein expression has been based on sequence optimisation that has favoured speed.
Nevertheless, pauses and slowdowns are pivotal to allow correct protein folding. 
Higher expression rates can be obtained by slowing down the translation process to allow a more efficient protein folding.
This tool has been developed to allow a protein sequence optimization that takes into account these aspects.

## Workflow
- Calculate codon usage: `python bin/dna2freq.py --cds data/organisms/ecoli/ncbi_dataset/data/GCF_000005845.2/cds_from_genomic.fna --organism ecoli`
- Equalise sequences: `ython bin/equalise.py --sequences test/seq.fasta`

## Data
**Note**: no exons/UTRs are present
- E. coli Genome assembly ASM584v2: cds sequences downloaded from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/
- H. sapiens assembly GRCh38.p14: cds sequences from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/


## TODO
- Finish equalise.py