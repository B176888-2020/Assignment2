#! /usr/bin/python3

################################## Modules and packages ##################################
# System modules
import os
import shutil
import subprocess
import sys

# String module
import string

# Regex
import re


################################## Functions ##################################
# `inputCheck` function is used to check and extract data from the inputs parameters
def inputCheck(proFamilyx, taxGroupx, spOrId = None, proSelection=None):
    # The function used to check the data type and extract data from inputs
    def extract_data(proFamilys):
        if type(proFamilys) == str:
            if re.search(".txt", proFamilys):
                fileProFamilys = open(proFamilys, "r")
                fileProFamilysContent = fileProFamilys.read()
                fileProFamilys.close()
                lsProFamilys = fileProFamilysContent.split("\n")
            else:
                lsProFamilys = proFamilys
        elif type(proFamilys) == list:
            lsProFamilys = proFamilys
        else:
            exit("\nError: The input format of protein families or taxonomy groups is wrong. \n" +
                 "The input of protein families should be a python list \n" +
                 "or string with pathway to the txt file containing protein familys \n" +
                 "or just the string of the protein family. \n" +
                 "Please check if the input format meets the requirement above.")
        return lsProFamilys

    # Extract the inputs and continue pairwise check
    lsProFamilys = extract_data(proFamilyx)
    lsTaxGroups = extract_data(taxGroupx)

    # Check if the datatype is pairwise
    if type(lsProFamilys) == type(lsTaxGroups):
        if type(lsProFamilys) == list:
            if len(lsProFamilys) != len(lsTaxGroups):
                exit("Error: The length of the protein families and taxonomy groups are not the same. \n" +
                     "Please check if your inputs have the pairwise protein families and taxonomy groups.")
    else:
        exit("Error: The data type of protein family and taxonomy group is different. \n " +
             "Please check if your inputs are both lists or strings")

    return [lsProFamilys, lsTaxGroups]


# `protree` function is used to analyse the conservation level of each protein family and taxonomy group
def protree(proFamily, taxGroup):
    # Get the primary protein sequence by `esearch` and `edirect`
    os.makedirs("./data", exist_ok=True)
    os.system("esearch -db protein -query '" +
              proFamily + "[protein] AND " + taxGroup + "[organism]' | efetch -db protein -format fasta > ./data/proSeq.fa")
    proSeq = open("./data/proSeq.fa", "r")
    proSeqContent = proSeq.read()
    proSeq.close()
    proSeqCount = proSeqContent.count(">")
    # TODO: print the result of different id/species
    # TODO: Summary the basic information of the sequences and report to the user
    # TODO: Checkpoint: inspect the data and check unusual data. Let user decide if they need to continue the process.
    # When the sequences is larger than the
    print("\nThere are " + str(proSeqCount) + " protein sequences of the protein family " + proFamily +
          "in the taxonomy group" + taxGroup + ". \n")

    if proSeqCount > 1000:
        proSeqCountReply = input("\nThe recommended number of sequences should be less than 1000, " +
                                 "otherwise the workload may influence the efficiency of this program. \n" +
                                 "Please enter YES or NO to decide whether you want to continue downstream analysis: \n")
        if proSeqCountReply == "NO":
            exit()

    # Choose the protein sequences by the sequence id or species names
    # Functions used for this section
    def species2id(speciesList, proSeqContent):
        speciesListr = speciesList.rstrip()
        pattern = '|'.join(speciesListr.split("\n"))
        sproID = open("./data/sproID.txt", "w")
        for sequence in proSeqContent.split(">"):
            if re.search(pattern, sequence):
                spId = sequence.split(" ")[0]
                sproID.write(spId + "\n")
            else:
                continue
        sproID.close()

    def id2seq(idList, modePrune):
        if re.search('.txt', idList):
            if re.search("EXCLUDE", modePrune):
                subprocess.call("pullseq -i ./data/proSeq.fa -n " + idList + " -e " + idList + " > ./data/proSeqS.fa",
                                shell=True)
            else:
                subprocess.call("pullseq -i ./data/proSeq.fa -n " + idList + " > ./data/proSeqS.fa", shell=True)
        else:
            exit("Error: The txt file containing the ID is not found.")

    # Select the sequences by species or protein id.
    proSeqCountReplyPrune = input(
        "\nDo you want to select target protein in specific species to conduct further analysis? \n" +
        "Please enter YES or NO to decide whether you want to select species: \n")
    if proSeqCountReplyPrune == "YES":
        modePrune = input(
            "\nIf you want to use protein sequence IDs to select protein sequences, please enter ID to search \n" +
            "If you want to use species name to select protein sequences, please enter SPECIES to search . \n" +
            "If you want to exclude these protein IDs or species names, you can add EXCLUDE before NORMAL or SPECIES symbol. \n" +
            "If you want to exit this prune section, you can enter EXIT to skip it. \n"
        )
        # subprocess.call("export PATH=/localdisk/data/BPSM/Assignment2/:$PATH")
        if re.search("ID", modePrune):
            idList = input(
                "\nPlease enter the protein IDs  that you want to select for further analysis and separate them with spaces . \n" +
                "If you have prepare the txt file containing protein id or species names separated by lines, you can provide the pathway of the file instead of typing them through the stdin: \n")
            if re.search('.txt', idList):
                id2seq(idList, modePrune)
            else:
                iproID = open("./data/iproID.txt", "w")
                lsid = idList.split(" ")
                for proid in lsid:
                    iproID.write(proid + "\n")
                iproID.close()
                id2seq("./data/iproID.txt", modePrune)
        elif re.search("SPECIES", modePrune):
            speciesList = input(
                "\nPlease enter the species names that you want to select for further analysis and separate them with \\n . \n" +
                "If you have prepare the txt file containing protein id or species names separated by lines, you can provide the pathway of the file instead of typing them through the stdin: \n")
            if re.search('.txt', speciesList):
                species = open(speciesList, "r")
                speciesContent = species.read()
                species.close()
                species2id(speciesContent, proSeqContent)
            else:
                speciesContent = '\n'.join(speciesList.split("\\n"))
                species2id(speciesContent, proSeqContent)
            id2seq("./data/sproID.txt", modePrune)
        elif re.search("EXIT", modePrune):
            print("The selection process has been skipped.")
        else:
            exit("Error: The command doesn't meet the requirement above. Please try again.")
    elif proSeqCountReplyPrune == "NO":
        subprocess.call("cp ./data/proSeq.fa ./data/proSeqN.fa", shell=True)

    # Align the protein sequences and determine the level of conservation by distance matrix
    print("Start Analysis")
    os.makedirs("./figures", exist_ok=True)
    os.system("clustalo --force --full --percent-id --maxnumseq 1000 --threads 12 -i ./data/proSeq?.fa " +
              "-o ./data/proAligned.fa --outfmt=fa --output-order=tree-order " +
              "--distmat-out=./data/proDistmat --guidetree-out=./data/proGuideTree")

    # Plot the guide tree
    os.system("figtree -graphic PDF ./data/proSeqGuideTree test.pdf")

    # Count the aligned sequence number in the file
    proAligned = open("./data/proAligned.fa", "r")
    proAlignedContent = proAligned.read()
    seqCount = proAlignedContent.count(">")
    proAligned.close()
    if seqCount > 250:
        seqCountReply = input("\nThere are " + str(seqCount) + " sequences in the sequence aligned file. \n" +
                              "It may influence further analysis if the number of sequences are more than 250. \n " +
                              "Please enter YES or NO to decide whether you want to continue the downstream analysis: \n")
        if seqCountReply == "NO":
            exit()
        elif seqCountReply == "YES":
            seqCountReply2 = input("Do you want to use the 250 most similar sequences as the input of plotcon? \n" +
                                   "Please enter YES or NO to decide whether you need this similarity prune: \n")
            if seqCountReply2 == "YES":
                proAlignedContent250 = '>'.join(proAlignedContent.split(">")[0:251])
                proAlignedPrune = open("./data/proAligned.fa", "w")
                proAlignedPrune.write(proAlignedContent250)
                proAlignedPrune.close()
            elif seqCountReply2 == "NO":
                exit()

    # Use `plotcon` program to visualise the conversation level
    os.system("plotcon -sformat fasta ./data/proAligned.fa -goutfile similarity -gdirectory ./figures")
    # Ask if user want to see the pictures
    eogReply = input("Do you want to see the visualisation results from plotcon directly? \n"
                     "Please enter YES or NO to decide whether you need to inspect the plotcon outputs: \n")
    if eogReply == "YES":
        os.system("eog ./figures/similarity*")

    # Scan protein sequence(s) of interest with motifs from the PROSITE database
    os.makedirs("./motifResults", exist_ok=True)
    os.system("patmatmotifs -sformat fasta -sequence ./data/proSeq.fa " +
              "-outfile motifReport -rdirectory2 ./motifResults")
    motifsReport = open("./motifResults/motifReport", "r")
    motifsReportContent = motifsReport.read()
    motifsReport.close()
    print("\n".join(element for element in motifsReportContent.split("\n") if 'Motif = ' in element))

    # Output the result


################################## Main program ##################################
# Get proper inputs from users
if len(sys.argv) == 1:
    reply = input("Protein family and taxonomic group are needed for proTree.py. \n" +
                  "If you just want to check the sample shown in the manual, please enter SAMPLE to continue. \n" +
                  "If you want to check your own datasets, please enter -v to follow an interactive process. \n" +
                  "If you want to exit, you can enter EXIT or press the Ctrl+C to exit the programme. \n")
    if reply == "EXIT":
        exit()
    elif reply == "SAMPLE":
        proFamilys = "glucose-6-phosphatase"
        taxGroups = "Aves"
    elif reply == "-v":
        proFamilys = input("Please enter the protein family: \n")
        taxGroups = input("Please enter the taxonomic group: \n")
elif sys.argv[1] == "-v":
    proFamilys = input("Please enter the protein family: \n")
    taxGroups = input("Please enter the taxonomic group: \n")
elif len(sys.argv) > 1:
    proFamilys = sys.argv[1]
    taxGroups = sys.argv[2]
    if len(sys.argv) > 3:
        spOrId = sys.argv[3]
        selectionls = sys.argv[4]  # absolute paths
    elif:
        spOrId = None
        selectionls = None

# Verify the input and Check whether the input values are suitable
inputCheckResult = inputCheck(proFamilys, taxGroups)
lsProFamilys = inputCheckResult[0]
lsTaxGroups = inputCheckResult[1]

# Conduct the conservation analysis
if type(lsProFamilys) == str:
    proFamily = lsProFamilys
    taxGroup = lsTaxGroups
    protree(proFamily, taxGroup)
else:
    for counter in list(range(0, len(lsProFamilys))):
        proFamily = lsProFamilys[counter]
        taxGroup = lsTaxGroups[counter]
        proSelection = selectionls[counter]
        dirPro = "./" + str(proFamily) + "_" + str(taxGroup)
        os.makedirs(dirPro, exist_ok=True)
        os.system("cd " + dirPro)
        protree(proFamily, taxGroup, spOrId, proSelection)
