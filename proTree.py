#! /usr/bin/python3

# Import modules and packages
# System modules
import os
import shutil
import subprocess
import sys
# String module
import string
# Regex
import re


# Get proper inputs from users
if len(sys.argv) == 1:
    reply = input("Protein family and taxonomic group are needed for proTree.py. \n" +
                  "If you just want to check the sample shown in the manual, please enter SAMPLE to continue. \n" +
                  "If you want to check your own datasets, please enter -v to follow an interactive process. \n" +
                  "If you want to exit, you can enter EXIT or press the Ctrl+C to exit the programme. \n")
    if reply == "EXIT":
        exit()
    elif reply == "SAMPLE":
        proFamily = "glucose-6-phosphatase"
        taxGroup = "Aves"
    elif reply == "-v":
        proFamily = input("Please enter the protein family: \n")
        taxGroup = input("Please enter the taxonomic group: \n")
elif sys.argv[1] == "-v":
    proFamily = input("Please enter the protein family: \n")
    taxGroup = input("Please enter the taxonomic group: \n")
elif len(sys.argv) > 1:
    proFamily = sys.argv[1]  # Generally, inputs from `sys.argv` will be settled as string
    taxGroup = sys.argv[2]  # The sys.argv[0] is defined as the name of the python script


# Get the primary protein sequence by `esearch` and `edirect`
os.makedirs("./data")
os.system("esearch -db protein -query '" +
          proFamily + "[protein] AND " + taxGroup + "[organism]' | efetch -db protein -format fasta > ./data/proSeq.fa")
# TODO: Summary the basic information of the sequences and report to the user
# TODO: Checkpoint: inspect the data and check unusual data. Let user decide if they need to continue the process.

# Determine, and plot the level of conservation
os.makedirs("./figures")
os.system("clustalo --force -i ./data/proSeq.fa -o ./data/proSeqAligned.fa --outfmt=fa --maxnumseq 1000 --threads 12")
os.system("plotcon -sformat fasta ./data/proSeqAligned.fa -graph png -goutfile similarity -gdirectory ./figures")
os.system("eog ./figures/similarity.1.png")

# Scan protein sequence(s) of interest with motifs from the PROSITE database
os.makedirs("./motifResults")
# TODO: check if prosextract is needed by locate function
os.system("patmatmotifs -sformat fasta -sequence ./data/proSeq.fa -outfile motifReport -rdirectory2 ./motifResults")


# Output the result


