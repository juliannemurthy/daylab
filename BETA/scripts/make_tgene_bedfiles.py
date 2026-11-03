""" Script that parses binding data from column 4 (RE locus) from a T-gene analysis file and separates it into three
columns: chrom, chromStart, chromEnd, stores it in a tab-delimited no header BED file for use in BETA. The T-gene
files are filtered by distance, keeping only peaks within 1 kb of a transcription start site (TSS). All the genes in
the files are significant by p-value too. The purpose is to pull out peaks we know are close to a TSS that may be buried
in the analysis, store into new BED files, and use those with BETA plus analysis.
 """
import pandas as pd
import argparse
import sys

def get_cli_args():
    """ Get command line options using argparse """
    # add parser description and define new flag -i
    parser = argparse.ArgumentParser(
        description="Provide T-gene analysis file to extract peak information into a BED file",
        )
    parser.add_argument('-i',
                        '--infile',
                        dest='INFILE',
                        type=str,
                        help='Path to the T-gene file to open',
                        required=True)
    parser.add_argument('-o',
                        '--outfile',
                        dest='OUTFILE',
                        type=str,
                        help='Path to the BED file to write to',
                        required=True)
    return parser.parse_args()

def parse_tgene_file(infile):
    """ Parses T-gene file for the fourth column 'RE locus', separates into three columns """
    # Load in the T-gene file
    tgene_df = pd.read_csv(infile, sep='\t')

    # Normalize column names strip spaces and make lowercase
    tgene_df.columns = tgene_df.columns.str.strip().str.lower()

    # Check if column of interest exists
    if "re_locus" not in tgene_df.columns:
        sys.exit("Error: Column 'RE_Locus' (or similar name) not found in file.")

    # Split RE_Locus into chrom, start, end. Splits if there is a : or -, and expands the split parts into separate DataFrame columns instead of Series of lists
    chrom_df = tgene_df["re_locus"].str.split("[:-]", expand=True)
    chrom_df.columns = ["chrom", "chromStart", "chromEnd"]

    return chrom_df

def write_bedfile(df, outfile):
    """ Stores chrom coordinate columns into a tab delimited header-less BED file """
    df.to_csv(outfile, sep='\t', index=False, header=False)


def main():
    """ Business Logic """
    # define args input variables
    args = get_cli_args()
    infile = args.INFILE
    outfile = args.OUTFILE
    # Parse the infile and separate the chromosome coordinate info
    chrom_df = parse_tgene_file(infile)
    # Write the info to BED file fit for BETA
    write_bedfile(chrom_df, outfile)

if __name__ == '__main__':
    main()
