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
os.makedirs("./data", exist_ok=True)
os.system("esearch -db protein -query '" +
          proFamily + "[protein] AND " + taxGroup + "[organism]' | efetch -db protein -format fasta > ./data/proSeq.fa")
# TODO: Summary the basic information of the sequences and report to the user
# TODO: Checkpoint: inspect the data and check unusual data. Let user decide if they need to continue the process.
os.system(proSeq.fa)
# Determine, and plot the level of conservation
os.makedirs("./figures", exist_ok=True)
os.system("clustalo --force --full --percent-id --maxnumseq 1000 --threads 12 -i ./data/proSeq.fa " +
          "-o ./data/proSeqAligned.fa --outfmt=fa --output-order=tree-order " +
          "--distmat-out=./data/proSeqDistmat --guidetree-out=./data/proSeqGuideTree")



# figtree -graphic PDF ./data/proSeqGuideTree test.pdf


# Count the aligned sequence number in the file
proSeqAligned = open("./data/proSeqAligned.fa", "r")
proSeqAlignedContent = proSeqAligned.read()
seqCount = proSeqAlignedContent.count(">")
proSeqAligned.close()
if seqCount > 250:
    seqCountReply = input("There are " + str(seqCount) + " sequences in the sequence aligned file. \n" +
                          "It may influence further analysis if the number of sequences are more than 250. \n " +
                          "Please enter YES or NO to decide whether you want to continue the downstream analysis: \n")
    if seqCountReply == "NO":
        exit()
    elif seqCountReply == "YES":
        seqCountReply2 = input("Do you want to use the 250 most similar sequences to plotcon? \n" +
                               "Please enter YES or NO to decide whether you need this similarity prune: \n")
        if seqCountReply2 == "YES":
            proSeqAlignedContent250 = '>'.join(proSeqAlignedContent.split(">")[0:251])
            proSeqAlignedPrune = open("./data/proSeqAligned.fa", "w")
            proSeqAlignedPrune.write(proSeqAlignedContent250)
            proSeqAlignedPrune.close()
        elif seqCountReply2 == "NO":
            exit()

os.system("plotcon -sformat fasta ./data/proSeqAligned.fa -goutfile similarity -gdirectory ./figures")
os.system("eog ./figures/similarity*")

# Scan protein sequence(s) of interest with motifs from the PROSITE database
os.makedirs("./motifResults", exist_ok=True)
# TODO: check if prosextract is needed by locate function
os.system("patmatmotifs -sformat fasta -sequence ./data/proSeq.fa " +
          "-outfile motifReport -rdirectory2 ./motifResults")
motifsReport = open("./motifResults/motifReport", "r")
motifsReportContent = motifsReport.read()
motifsReport.close()

print()


# Output the result


