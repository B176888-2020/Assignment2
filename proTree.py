#! /usr/bin/python3

################################## Modules and packages ##################################
# System modules
import os
import subprocess
import sys

# String
import string

# Regex
import re

################################## Manual ##################################
manual = "\nNAME\n" \
         "    proTree: analyse the conservation level of protein family within taxonomy group\n" \
         "SYNOPSIS\n" \
         "    proTree [-options] [protein_families] [taxonomy_groups] [project_directory] [selection_modes] [selection]\n" \
         "DESCRIPTION\n" \
         "    proTree can analyse and visualise the conservation level of the protein family within taxonomy group.\n" \
         "    Functions in EMBOSS are also provdied to broaden the scale of motifs search and other usages.\n" \
         "\n" \
         "    Options:\n" \
         "\n" \
         "    -h, --help\n" \
         "        display this help and exit\n" \
         "    -v\n" \
         "        use interactive mode to conduct the analysis process. Users can provide inputs acccording to the requirements.\n" \
         "    -s\n" \
         "        use silent mode to conduct the analysis process if users can provide correct inputs.\n" \
         "\n" \
         "AUTHOR\n" \
         "    B176888-2020\n" \
         "REPORTINGS BUGS\n" \
         "    Please feel free to Report bugs to <https://github.com/B176888-2020/Assignment2/issues>\n" \
         "SEE ALSO\n" \
         "    Full program and documentation at <https://github.com/B176888-2020/Assignment2>."


################################## Functions ##################################
# `inputCheck` function is used to check and extract data from the inputs parameters
def inputCheck(proFamilyx, taxGroupx, spOrId=None, proSelection=None):
    # Extract the data from differet types of inputs
    def extract_data(proFamilys):
        if re.search('\[', proFamilys):
            proFamilys = proFamilys.strip('[]').split(',')
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

    # Extract the inputs and continue the pairwise check
    lsProFamilys = extract_data(proFamilyx)
    lsTaxGroups = extract_data(taxGroupx)
    if spOrId is not None:
        lsSpOrId = extract_data(spOrId)
    else:
        lsSpOrId = None
    if proSelection is not None:
        lsproSelection = extract_data(proSelection)
    else:
        lsproSelection = None

    # Check if the datatype is pairwise
    if type(lsProFamilys) == type(lsTaxGroups):
        if type(lsProFamilys) == list:
            if len(lsProFamilys) != len(lsTaxGroups):
                exit("\nError: The length of the protein families and taxonomy groups are not the same. \n" +
                     "Please check if your inputs have the pairwise protein families and taxonomy groups.")
    else:
        exit("\nError: The data type of protein family and taxonomy group is different. \n " +
             "Please check if your inputs are both lists or strings")

    return lsProFamilys, lsTaxGroups, lsSpOrId, lsproSelection


# `protree` function is used to analyse the conservation level of each protein family and taxonomy group
def protree(proFamily, taxGroup, spOrId, proSelection, vbo):
    # Mark the start of the analysis
    print("\n################################## Data Collection and Selection ##################################")

    # Get the primary protein sequence by `esearch` and `edirect`
    print("\n**** Collecting the data by esearch and efetch... ****")
    os.makedirs("./data", exist_ok=True)
    os.makedirs("./sum_data", exist_ok=True)
    os.system("esearch -db protein -query '" +
              proFamily + " AND " + taxGroup +
              "[organism] NOT Partial' | efetch -db protein -format fasta > ./data/proSeq.fa")
    print("Done\n")

    # Get the sequence quantity and the number of species
    proSeq = open("./data/proSeq.fa", "r")
    proSeqContent = proSeq.read()
    proSeq.close()
    proSeqCount = proSeqContent.count(">")
    lsp = []
    for sequence in proSeqContent.split(">")[1:]:
        if re.search("\[(.*)\]", sequence):
            spn = re.search("\[(.*)\]", sequence).group(1)
            lsp.append(spn)
        else:
            continue
    lsp = set(lsp)
    # Print the numbers of id/species in the dataset
    print("\nBasic Information of collected sequences:")
    print("Sequences Quantity: " + str(proSeqCount) + "\n" +
          "Species Quantity: " + str(len(lsp)) + "\n")
    if len(lsp) < 10:
        print("Species: \n")
        for sp in lsp:
            print(str(sp))
    else:
        print("The number of species is more than 10. Thus the names of these species will not be shown on stdout.\n")

    # Use pepstats to calculate statistics of protein sequecnes' properties
    print("\n**** Summarise the protein properties by pepstats method... ****")
    if vbo:
        pepReply = input(
            "\nDo you want to calculate statistics of protein sequecnes' properties by pepstats method? \n" +
            "Please enter YES or NO to decide whether you want to use pepstats method: \n")
        if pepReply.upper() == "YES":
            os.system("pepstats -sequence ./data/proSeq.fa --outfile ./sum_data/proSeqProperties.txt")
            print("Done. The output will be stored in the proSeqProperties.txt in the sum_data directory. \n")
        else:
            print("Skip: pepstats")
    else:
        os.system("pepstats -sequence ./data/proSeq.fa --outfile ./sum_data/proSeqProperties.txt")
        print("Done. The output will be stored in the proSeqProperties.txt in the sum_data directory. \n")

    # Give warnings when the number of sequences is more than 1000
    if proSeqCount > 1000:
        if vbo:
            proSeqCountReply = input("\nThe recommended number of sequences should be less than 1000, \n" +
                                     "otherwise it may take a long time to run through this program. \n" +
                                     "The selection methods provided in the following part may help you select proper species or id to reduce the workload \n"
                                     "Please enter YES or NO to decide whether you want to continue downstream analysis: \n")
            if proSeqCountReply.upper() == "NO":
                exit()
        else:
            print("\nWarning: The number of total sequences is larger than 1000.\n")

    # Choose the protein sequences by the sequence id or species names
    print("\n**** Select sequences by protein sequence id or species names... ****")

    # Transfer species to protein sequence id
    def species2id(speciesList, proSeqContent):
        speciesListr = speciesList.rstrip()
        pattern = '|'.join(speciesListr.split("\n"))
        sproID = open("./data/sproID.txt", "w")
        for sequence in proSeqContent.split(">")[1:]:
            if re.search(pattern, sequence):
                spId = sequence.split(" ")[0]
                sproID.write(spId + "\n")
            else:
                continue
        sproID.close()

    # Extract the protein sequences according to protein sequence id
    def id2seq(idList, modePrune):
        if re.search('.txt', idList):
            if re.search("EXCLUDE", modePrune):
                subprocess.call("pullseq -i ./data/proSeq.fa -n " + idList + " -e " + idList + " > ./data/proSeqS.fa",
                                shell=True)
            else:
                subprocess.call("pullseq -i ./data/proSeq.fa -n " + idList + " > ./data/proSeqS.fa", shell=True)
        else:
            exit("\nError: The txt file containing the ID is not found.")

    # The function to transfer and select protein sequence id and extract corresponding sequences
    def seqSelection(ls, modePrune):
        if re.search("ID", modePrune):
            if re.search('.txt', ls):
                id2seq(ls, modePrune)
            else:
                iproID = open("./data/iproID.txt", "w")
                lsid = ls.split('\\t')
                for proid in lsid:
                    iproID.write(proid + "\n")
                iproID.close()
                id2seq("./data/iproID.txt", modePrune)
        elif re.search("SP", modePrune):
            if re.search('.txt', ls):
                species = open(ls, "r")
                speciesContent = species.read()
                species.close()
                species2id(speciesContent, proSeqContent)
            else:
                speciesContent = '\n'.join(ls.split("\\n"))
                species2id(speciesContent, proSeqContent)
            id2seq("./data/sproID.txt", modePrune)
        else:
            print("Error: The selection mode and document may be not matched.")

    # Select the sequences by species or protein id in different mod
    if (spOrId is None and proSelection is not None) or (proSelection is None and spOrId is not None):
        exit("\nError: One of proSelection or spOrId is not found. These two inputs should be both provided or provided neither.")
    elif (spOrId is None) and (proSelection is None):
        if vbo:
            proSeqCountReplyPrune = input(
                "\nDo you want to select target protein in specific species to conduct further analysis? \n" +
                "Please enter YES or NO to decide whether you want to select species: \n")
            if proSeqCountReplyPrune.upper() == "YES":
                modePrune = input(
                    "\nIf you want to use protein sequence IDs to select protein sequences, please enter ID\n" +
                    "If you want to use species name to select protein sequences, please enter SP. \n" +
                    "If you want to exclude these protein IDs or species names, you can add EXCLUDE before ID or SP symbol. \n" +
                    "If you want to skip this selection process, you can enter SKIP to skip it. \n"
                )
                modePrune = modePrune.upper()

                if re.search("ID", modePrune):
                    idList = input(
                        "\nPlease enter the protein IDs that you want to select for further analysis and separate them with \\t . \n" +
                        "If you have prepare the txt file containing protein id or species names separated by lines, you can also provide the pathway to the file: \n")
                    seqSelection(idList, modePrune)
                    selected = True
                    print("Done\n")
                elif re.search("SP", modePrune):
                    speciesList = input(
                        "\nPlease enter the species names that you want to select for further analysis and separate them with \\n . \n" +
                        "If you have prepare the txt file containing protein id or species names separated by lines, you can also provide the pathway to the file: \n")
                    seqSelection(speciesList, modePrune)
                    selected = True
                    print("Done\n")
                elif re.search("SKIP", modePrune):
                    print("Skip: The selection process has been skipped.")
                    subprocess.call("cp ./data/proSeq.fa ./data/proSeqN.fa", shell=True)
                    selected = False
                else:
                    exit("\nError: The command doesn't meet the requirement above. Please try again.")
            elif proSeqCountReplyPrune.upper() == "NO":
                print("Skip: The selection process has been skipped.")
                subprocess.call("cp ./data/proSeq.fa ./data/proSeqN.fa", shell=True)
                selected = False
        else:
            subprocess.call("cp ./data/proSeq.fa ./data/proSeqN.fa", shell=True)
            selected = False
    elif (spOrId is not None) and (proSelection is not None):
        modePrune = spOrId.upper()
        print("\nSelect sequences according to " + str(proSelection) + "...")
        seqSelection(proSelection, modePrune)
        selected = True
        print("Done\n")
    else:
        exit(
            "\nError: The options for species or protein id to select sequences may not meet the requirements. \n" +
            " Please check your inputs and try again.")

    # Count the sequnces after the selection
    if selected:
        proHth = open("./data/proSeqS.fa", "r")
    else:
        proHth = open("./data/proSeqN.fa", "r")
    proHthContent = proHth.read()
    proHth.close()
    proHthCount = proHthContent.count(">")
    lsphth = []
    for sequence in proHthContent.split(">")[1:]:
        spHth = re.search('\[(.*)\]', sequence).group(1)
        lsphth.append(spHth)
    lsphth = set(lsphth)
    # Print the result of different id/species
    print("\nBasic Information of chosen sequences:")
    print("Sequences Quantity: " + str(proHthCount) + "\n" +
          "Species Quantity: " + str(len(lsphth)) + "\n")
    if len(lsphth) < 10:
        print("Species:")
        for sp in lsphth:
            print(str(sp))
    else:
        print("The number of species is more than 10. Thus the names of these species will not be shown on stdout.\n")

    # Align the protein sequences and determine the level of conservation by distance matrix
    print("\n################################## MSA and Conservation ##################################")

    # Check if the data set is too small to do MSA
    if proHthCount < 2:
        exit("Error: The number of protein sequences is too small to conduct multiple sequence alignment.")
    # MSA process and generate distance matrix
    print("\n**** Processing multiple sequence alignment... ****")
    os.makedirs("./figures", exist_ok=True)
    os.system("clustalo --force --full --percent-id --threads 12 -i ./data/proSeq?.fa " +
              "-o ./data/proAligned.fa --outfmt=fa --output-order=tree-order " +
              "--distmat-out=./data/proDistmat --guidetree-out=./data/proGuideTree")
    print("Done\n")

    # Summarise the aligned data information
    print("\n**** Using infoalign to summarise the protein sequences... ****")
    if vbo:
        infoalign_reply = input("\nDo you want to get the summary of the protein sequences by infoalign? \n" +
                                "Please enter YES or NO to decide whether you want to use infoalign method: \n")
        if infoalign_reply.upper() == "YES":
            os.system("infoalign -sequence ./data/proAligned.fa -outfile ./sum_data/proInfoAlign.infoalign")
            print("Done\n")
        elif infoalign_reply.upper() == "NO":
            print("\nSkip: infoalign")
    else:
        os.system("infoalign -sequence ./data/proAligned.fa -outfile ./sum_data/proInfoAlign.infoalign")
        print("Done\n")


    # Plot the guide tree
    print("\n**** Using figtree to plot the guide tree... ****")
    if vbo:
        guidetree_reply = input("\nDo you want to plot the guide tree by figtree? \n" +
                                "Please enter YES or NO to decide whether you want to use figtree method: \n")
        if guidetree_reply.upper() == "YES":
            os.system("figtree -graphic PDF ./data/proGuideTree ./figures/guide_tree.pdf")
            print("Done\n")
        elif guidetree_reply.upper() == "NO":
            print("\nSkip: figtree")
    else:
        os.system("figtree -graphic PDF ./data/proGuideTree ./figures/guide_tree.pdf")
        print("Done\n")

    # Test if the quantity of the sequences is larger than 250 and whether to prune the sequences by similarity
    def seq250(proAlignedContent):
        proAlignedContent250 = '>'.join(proAlignedContent.split(">")[0:251])
        proAlignedPrune = open("./data/proAligned.fa", "w")
        proAlignedPrune.write(proAlignedContent250)
        proAlignedPrune.close()

    print("\n**** Check if quantity of the sequences is too large for plotcon... ****")
    proAligned = open("./data/proAligned.fa", "r")
    proAlignedContent = proAligned.read()
    seqCount = proAlignedContent.count(">")
    proAligned.close()
    if vbo:
        if seqCount > 250:
            seqCountReply = input("\nThere are " + str(seqCount) + " sequences in the sequence aligned file. \n" +
                                  "It may influence the plotcon function if the number of sequences are more than 250. \n" +
                                  "Please enter YES or NO to decide whether you want to continue the downstream analysis: \n")
            if seqCountReply.upper() == "NO":
                exit()
            elif seqCountReply.upper() == "YES":
                seqCountReply2 = input(
                    "\nDo you want to use the 250 most similar sequences as the input of plotcon? \n" +
                    "Please enter YES or NO to decide whether you need this similarity prune: \n")
                if seqCountReply2.upper() == "YES":
                    print("\n**** Select the 250 most similar sequences as the input of plotcon... ****")
                    seq250(proAlignedContent)
                    print("Done\n")
                elif seqCountReply2.upper() == "NO":
                    exit()
    else:
        print("\n**** Select the 250 most similar sequences as the input of plotcon... ****")
        seq250(proAlignedContent)
        print("Done\n")

    # Use `plotcon` program to visualise the conversation level
    print("\n**** Using plotcon to visualise the conversation level of protein sequences... ****")
    if vbo:
        os.system(
            "plotcon -sformat fasta ./data/proAligned.fa -sprotein1 Yes -goutfile similarity -gdirectory ./figures")
        eogReply = input("Do you want to see the visualisation result from plotcon directly? \n"
                         "Please enter YES or NO to decide whether you need to inspect the plotcon outputs: \n")
        if eogReply.upper() == "YES":
            os.system("eog ./figures/similarity*")
    else:
        os.system(
            "plotcon -sformat fasta ./data/proAligned.fa -sprotein1 Yes -goutfile similarity -gdirectory ./figures -winsize 4 -graph svg")
    print("Done\n")

    # Search different types of motifs
    print("\n################################## Motifs and Special Sites ##################################")
    # Scan protein sequence(s) of interest with motifs from the PROSITE database
    print("\n**** Searching the motifs from the PROSITE database... ****")
    os.makedirs("./motifs", exist_ok=True)
    os.makedirs("./motifResults", exist_ok=True)
    os.makedirs("./motifResults2", exist_ok=True)
    if vbo:
        motifsReply = input("\nDo you want to show the motifs information in stdout? \n" +
                            "Please enter YES or NO to choose whether to show the motifs.\n")
    else:
        motifsReply = "NO"
    if motifsReply.upper() == "YES":
        print("\nMotifs: \n")
    for seq in proHthContent.split(">")[1:]:
        seqId = seq.split(" ")[0]
        if motifsReply.upper() == "YES":
            print("\n** " + seqId + " **")
        seqfilename = "./motifs/" + seqId + ".fa"
        seqfile = open(seqfilename, "w")
        seqfile.write(">" + seq)
        seqfile.close()
        if motifsReply.upper() == "YES":
            os.system("patmatmotifs -sformat fasta -sequence " + seqfilename +
                      " -outfile " + seqId + "_motifs_report.patmatmotifs -rdirectory2 ./motifResults")
        else:
            os.system("patmatmotifs -auto Y -sformat fasta -sequence " + seqfilename +
                      " -outfile " + seqId + "_motifs_report.patmatmotifs -rdirectory2 ./motifResults")
        if motifsReply.upper() == "YES":
            seqmotiffilename = "./motifResults/" + seqId + "_motifs_report.patmatmotifs"
            motifsReport = open(seqmotiffilename, "r")
            motifsReportContent = motifsReport.read()
            motifsReport.close()
            motifs = "\n".join(element for element in motifsReportContent.split("\n") if 'Motif = ' in element)
            print(motifs)
    print("Done. The output results will be stored in ./motifResults/motifsReport\n")

    # Additional motifs function discovery:
    # Find potentially antigenic regions of protein sequences
    print("\n**** Find potentially antigenic regions of a protein sequence... ****")
    if vbo:
        antiReply = input("\nDo you want to find potentially antigenic regions of a protein sequence by antigenic? \n" +
                          "Please enter YES or NO to decide whether you need to use antigenic program: \n")
        if antiReply.upper() == "YES":
            os.system("antigenic -sequence ./data/proSeq?.fa -minlen 6 -outfile ./motifResults2/antigenic_sites.antigenic")
            print("Done. The output will been stored in ./motifResults2/antigenic_sites.antigenic\n")
        elif antiReply.upper() == "NO":
            print("Skip: antigenic")
    else:
        os.system("antigenic -sequence ./data/proSeq?.fa -minlen 6 -outfile ./motifResults2/antigenic_sites.antigenic")
        print("Done. The output will been stored in ./motifResults2/antigenic_sites.antigenic\n")

    # Find helix-turn-helix nucleic acid binding motifs
    def hth(proHthContent):
        proHthContentr = proHthContent.replace("BZ", "-").replace("U", "-").replace("X", "-")
        proHth2 = open("./data/proHth.fa", "w")
        proHth2.write(proHthContentr)
        proHth2.close()
        os.system("helixturnhelix -sequence ./data/proHth.fa -outfile ./motifResults2/n2p.hth")

    print("\n**** Find helix-turn-helix nucleic acid binding motifs... ****")
    if vbo:
        hthReply = input("\nDo you want to find find helix-turn-helix nucleic acid binding motifs by helixturnhelix? \n" +
                          "Please enter YES or NO to decide whether you need to use helixturnhelix program: \n")
        if hthReply.upper() == "YES":
            hth(proHthContent)
            print("Done.  The output will been stored in ./motifResults2/n2p.hth\n")
        elif hthReply.upper() == "NO":
            print("Skip: helixturnhelix")
    else:
        hth(proHthContent)
        print("Done.  The output will been stored in ./motifResults2/n2p.hth\n")
    print("\n The Data Process for " + proFamily + " and " + taxGroup + ": Finish.")

# The major analysis process for every pair/group of inputs
def main(proFamily, taxGroup, projectSpace, lsSpOrId, lsproSelection, vbo):
    dirPro = projectSpace + str(proFamily.replace(" ", "-")) + "_" + str(taxGroup.replace(" ", "-"))
    os.makedirs(dirPro, exist_ok=True)
    os.chdir(dirPro)
    print("\nInput information: \n" +
          "Protein Family: " + str(proFamily) + "\n" +
          "Taxonomy Group: " + str(taxGroup) + "\n" +
          "Selection Mode: " + str(lsSpOrId) + "\n" +
          "Selection Target: " + str(lsproSelection) + "\n" +
          "Project Location: " + str(projectSpace) + "\n" +
          "Output Directory: " + str(dirPro)
          )
    protree(proFamily, taxGroup, lsSpOrId, lsproSelection, vbo)
    os.chdir("../")


################################## Main program ##################################

# server: PATH for esearch and efetch
os.system("echo \"\nexport PATH=$PATH:/localdisk/home/$USER/edirect/\" >> $HOME/.bashrc")
os.system("echo \"\nexport PATH=$PATH:/localdisk/data/BPSM/Assignment2/\" >> $HOME/.bashrc")
os.system("source $HOME/.bashrc")

# Some default values for the input arguments
projectSpace = "./"
spOrId = None
selectionls = None
vbo = True

# The start marker
print("\n################################## proTree ##################################")

# Get proper inputs from users
if len(sys.argv) == 1:
    reply = input("\nInteractive Mode(Default): \n" +
                  "If you want to check the sample shown in the manual, please enter SAMPLE to continue. \n" +
                  "If you want to check your own datasets, please enter -v to follow an interactive process. \n" +
                  "More detailed information about proTree usage can be checked by enter -h or --help. \n" +
                  "If you want to exit, please enter EXIT or press the Ctrl+C to exit the programme. \n")
    if reply == "EXIT":
        exit()
    elif reply == "-h" or reply == "--help":
        print(manual)
        exit()
    elif reply == "SAMPLE":
        proFamilys = "glucose-6-phosphatase"
        taxGroups = "Aves"
    elif reply == "-v":
        proFamilys = input("\nPlease enter the protein family: \n")
        taxGroups = input("\nPlease enter the taxonomic group: \n")
elif len(sys.argv) > 1:
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(manual)
        exit()
    elif sys.argv[1] == "-v" and len(sys.argv) == 2:
        print("\nInteractive Mode(Default):")
        proFamilys = input("\nPlease enter the protein family: \n")
        taxGroups = input("\nPlease enter the taxonomic group: \n")
    elif sys.argv[1] == "-v" or sys.argv[1] == "-s":
        if sys.argv[1] == "-v":
            print("\nInteractive Mode(Default):")
        elif sys.argv[1] == "-s":
            vbo = False
            print("\nSilent Mode:")
        proFamilys = sys.argv[2]
        taxGroups = sys.argv[3]
        if len(sys.argv) > 4:
            projectSpace = sys.argv[4]
        if len(sys.argv) > 5:
            spOrId = sys.argv[5]
        if len(sys.argv) > 6:
            selectionls = sys.argv[6]  # absolute paths
    else:
        print("\nInteractive Mode(Default): \n")
        proFamilys = sys.argv[1]
        taxGroups = sys.argv[2]
        if len(sys.argv) > 3:
            projectSpace = sys.argv[3]
        if len(sys.argv) > 4:
            spOrId = sys.argv[4]
        if len(sys.argv) > 5:
            selectionls = sys.argv[5]  # absolute paths

# Verify the input and Check whether the input values are suitable
print("\n**** Validating the inputs... ****")
lsProFamilys, lsTaxGroups, lsSpOrId, lsproSelection = inputCheck(proFamilys, taxGroups, spOrId, selectionls)
print("Done")

# Conduct the conservation analysis
if type(lsProFamilys) == str:
    proFamily = lsProFamilys
    taxGroup = lsTaxGroups
    main(proFamily, taxGroup, projectSpace, lsSpOrId, lsproSelection, vbo)
elif type(lsProFamilys) == list:
    for counter in list(range(0, len(lsProFamilys))):
        print("\n################################## Total Process: " + str(counter + 1) + "/" + str(len(lsProFamilys)) + " ##################################\n")
        proFamily = lsProFamilys[counter]
        taxGroup = lsTaxGroups[counter]
        if lsSpOrId is not None:
            if type(lsSpOrId) == str:
                sporid = lsSpOrId
            else:
                sporid = lsSpOrId[counter]
        else:
            sporid = None
        if lsproSelection is not None:
            if type(lsproSelection) == str:
                proSelection = lsproSelection
            else:
                proSelection = lsproSelection[counter]
        else:
            proSelection = None
        main(proFamily, taxGroup, projectSpace, sporid, proSelection, vbo)
