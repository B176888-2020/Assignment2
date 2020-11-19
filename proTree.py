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
                  "If you just want to check the default sample shown in the manual, please enter SAMPLE to continue. \n" +
                  "If you want to check your own datasets, please enter -i to follow an interactive process. \n" +
                  "If you want to exit, you can enter EXIT or press the Ctrl+C to exit the programme. \n")
    if reply == "EXIT":
        exit()
    elif reply == "SAMPLE":
        proFamily = "glucose-6-phosphatase"
        taxGroup = "Aves"
    elif reply == "-i":
        proFamily = input("Please enter the protein family. \n")
        taxGroup = input("Please enter the taxonomic group. \n")
elif sys.argv[1] == "-i":
    proFamily = input("Please enter the protein family. \n")
    taxGroup = input("Please enter the taxonomic group. \n")
elif len(sys.argv) > 1:
    proFamily = sys.argv[1]  # Generally, inputs from `sys.argv` will be settled as string
    taxGroup = sys.argv[2]  # The sys.argv[0] is defined as the name of the python script


# Testing variables
var1 = proFamily
var2 = taxGroup

# Output the result
print(var1 + var2)

