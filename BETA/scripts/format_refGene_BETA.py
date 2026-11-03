""" Script to format refGene annotation file for BETA analysis from UCSC refGene.txt file """
import pandas as pd
import gzip

# Load and parse the refGene annotation file
refgene_file_path = "/Users/juliannemurthy/refGene.txt.gz"
with gzip.open(refgene_file_path, 'rt') as file:
    refgene_df = pd.read_csv(file, sep="\t", header=None)

# Assign column names based on UCSC refGene.txt format
refgene_df.columns = [
    "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
    "exonCount", "exonStarts", "exonEnds", "score", "geneName",
    "cdsStartStat", "cdsEndStat", "exonFrames"
]

# Select required columns for BETA
refgene_beta = refgene_df[["name", "chrom", "strand", "txStart", "txEnd", "geneName"]]

# Drop NA missing values in any row
refgene_beta = refgene_beta.dropna()

# Sort by chromosome and start position
refgene_beta = refgene_beta.sort_values(by=["chrom", "txStart"])

# Save without header and tab-delimited
output_path = "/Users/juliannemurthy/refgene_for_BETA.bed"
refgene_beta.to_csv(output_path, sep='\t', index=False, header=False)


