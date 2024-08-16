#!/usr/bin/env python

"""
This module provides functions to extract sequences from a FASTA or FASTQ file based on a list of sequence IDs.
The sequence header in the input file must be in Uniprot format, i.e., ">sp|ID|GENE_NAME Protein Name".

Authors: Sebastian Gritsch, Karin Binder

It includes the following functions:

- check_format: Checks if a file is in FASTA or FASTQ format.
- parse_query: Parses a query file and returns a set of sequence IDs.
- get_subset_fast: Prints a subset of reads from a FASTA or FASTQ file based on IDs from a query file. This version is faster as it calls the correct parser directly.
- get_subset_parse: Prints a subset of reads from a FASTA or FASTQ file based on IDs from a query file. Note that this function does not print quality scores for "fastq" files.
- parse_arguments: Parses command line arguments.
- main: Main function that orchestrates the sub-sampling process.

The module uses the Biopython library for handling sequence files and includes logging to provide information during execution.

Usage:
$    extract_sequences.py -f FILENAME -q QUERY [-f] [-p]

Options:
  -h, --help            Show this help message and exit.
  -f FILENAME, --filename=FILENAME  Input file in FASTA or FASTQ format.
  -q QUERY, --query=QUERY  File with IDs to extract.
  -f, --fast            Use fast mode.
"""

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import sys
import logging
import time
import os
import pandas as pd

# Create and configure logger
# https://docs.python.org/3/howto/logging.html,
# https://stackoverflow.com/a/56144390/
logging.basicConfig(
    level=logging.NOTSET,
    filename="subsample.log",
    filemode="w",
    format="%(levelname)s:%(asctime)s:%(message)s",
)  # configure root logger
logger = logging.getLogger(__name__)  # create custom logger
# Logging levels: DEBUG/INFO/WARNING/ERROR/CRITICAL
logger.setLevel(logging.INFO)  # set logging level for our logger


# check file format
def check_format(filename):
    """Check if file is fasta or fastq and return the file format as string."""
    # check first line for @ or >
    with open(filename, "r") as f:
        line = f.readline()
        if line.startswith("@"):
            format = "fastq"
        elif line.startswith(">"):
            format = "fasta"
        else:
            error_message = f"{filename} is not a valid FASTA or FASTQ file"
            logger.error(error_message)
            raise ValueError(error_message)

    logger.info(f"File format: {format}")
    return format


# parse query file
def parse_query(query_file: str) -> set[str]:
    """Parse query file and return a list of sequence IDs."""
    df = pd.read_csv(query_file)
    query = set(df["Protein.ID"])

    logger.info(f"Number of IDs in query file: {len(query)}")
    return query


def get_subset_fast(filename, query, format):
    """
    Print a subset of reads from a FASTA or FASTQ file based on IDs from a query file. This version is faster as it calls the correct parser directly.

    Args:
        filename (str): The path to the input file.
        query (list): A list of sequence IDs to extract.
        format (str): The format of the input file ('fasta' or 'fastq').

    Returns:
        None
    """
    # Open file and iterate through records to get the total number of reads
    total_reads = 0
    with open(filename, "r") as file:
        if format == "fasta":
            for record in SimpleFastaParser(file):
                total_reads += 1
        elif format == "fastq":
            for record in FastqGeneralIterator(file):
                total_reads += 1

        logger.info(f"Total number of reads in {filename}: {total_reads}")

        # Jump to beginning of file
        file.seek(0)

        # Iterate through records again and print the selected ones
        if format == "fasta":
            while query:
                for record in SimpleFastaParser(file):
                    if record[0].startswith("rev"):  # skip reverse reads
                        continue
                    elif record[0].split("|")[1] in query:
                        print(f">{record[0]}")  # print sequence id
                        print(record[1], end="\n")  # print sequence
                        query.remove(record[0].split("|")[1])  # remove seqid from query
        elif format == "fastq":
            while query:
                for record in SimpleFastaParser(file):
                    if record[0].startswith("rev"):  # skip reverse reads
                        continue
                    elif record[0].split("|")[1] in query:
                        print(f"@{record[0]}")  # print sequence id
                        print(record[1])  # print sequence
                        print("+")  # print placeholder for decription line
                        print(record[2], end="\n")  # print quality score
                        query.remove(record[0].split("|")[1])  # remove seqid from query


def get_subset_parse(filename, query, format):
    """
    Print a subset of reads from a FASTA or FASTQ file based on IDs from a query file.
    This function only prints the sequence id followed by the sequence as a single string.

    Args:
        filename (str): The path to the input file.
        query (list): A list of sequence IDs to extract.
        format (str): The format of the input file. Must be either "fasta" or "fastq".

    Returns:
        None

    """
    #
    if format == "fastq":
        logger.warning(
            "SeqIO.parse cannot print quality scores. Use -f or -p option instead."
        )

    # Iterate through records and get total number of reads
    records = SeqIO.parse(filename, format)
    total_reads = sum(1 for _ in records)
    logger.info(f"Total number of reads in {filename}: {total_reads}")

    # Get record subsample
    while query:
        for record in SeqIO.parse(filename, format):
            if record.id.startswith("rev"):  # skip reverse reads
                continue
            elif record.id.split("|")[1] in query:
                if format == "fasta":
                    print(f">{record.description}")  # print sequence id
                    print(record.seq, end="\n")  # print sequence
                    query.remove(record.id.split("|")[1])  # remove seqid from query
                elif format == "fastq":
                    print(f"@{record.description}")  # print sequence id
                    print(record.seq, end="\n")  # print sequence
                    query.remove(record.id.split("|")[1])  # remove seqid from query


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract sequences from  aFASTA or FASTQ file based on a list of IDs."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="input file in FASTA or FASTQ format."
    )
    parser.add_argument(
        "-q", "--query", required=True, help="file with IDs to extract."
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "--fast", action="store_true", help="use fast mode")

    arguments = parser.parse_args()
    return arguments


def main(args):
    """main function"""
    start_time = time.time()
    logger.info(vars(args))  # log command line arguments
    filename = args.input

    # Check file format
    format = check_format(filename)

    # Parse query file
    logger.info(f"Parsing query file: {args.query}")
    query = parse_query(args.query)
    logger.info(f"Number of IDs in query file: {len(query)}")

    # Get subset
    logger.info(f"Getting subset of {len(query)} reads.")
    if args.fast:
        get_subset_fast(filename, query, format)
    else:
        get_subset_parse(filename, query, format)

    end_time = time.time()
    logger.info(f"Job done in {end_time - start_time:.2f} seconds.")


if __name__ == "__main__":
    arguments = parse_arguments()
    main(arguments)
