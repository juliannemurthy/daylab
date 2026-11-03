#!/usr/bin/env python


# run this Python code (in a Python interpreter) from a folder where all files are present
from bioinfokit.analys import HtsAna
# make sure all individual count files are present in same folder
# by default, it assumes each count file has .txt extension 
HtsAna.merge_featureCount()
