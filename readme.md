# Automated ordering of enzyme designs 

This repo is a collection of scripts that can be used together to automate
ordering of enzyme designs via Kunkel mutagenesis on Transcriptic. 

## How it works 

## Instructions for use 

Put your designed PDBs into a folder and modify `order.bash` to look for
designs in that folder. Also include a path to the wild type nucleotide
sequence of the gene that was used to generate the ssDNA. Run

```bash
bash order.bash | transcriptic submit -p project_name
```

to submit the Kunkel run. 

## How it works 

The script will execute a glob query to find all the designs in a folder
and iterate over them. It will pull out the protein sequence from the PDB
and align the protein sequence to the nucleotide sequence using `blastx`.
Python helper scripts then diff the two sequences, generate the mutagenic
oligos necessary to transform one sequence into another, format the oligos
for ordering, and generate autoprotocol for use on Transcriptic's platform. 

Questions, comments, and pull requests welcome as this is still very much
a work in progress! 
