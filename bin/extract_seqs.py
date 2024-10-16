#!/usr/bin/env python3

##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
#

"""Provide a command line tool to fetch a list of refseq genome ids to a single file, useful for kraken2 database building or alignment purposes"""

import argparse
import sys
import os

from Bio import SeqIO
import gzip

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Import a list of taxids, find all classified k2 reads from outfile with at least one occurrence of that taxid",
    )
    parser.add_argument(
        "-i",
        "--taxids",
        required=True,
        metavar="TAXIDS",
        help="Taxid list, one taxid per line",
    )
    parser.add_argument(
        "-k",
        "--k2",
        required=True,
        metavar="K2OUT",
        help="Kraken2 output per-read file",
    )
    parser.add_argument(
        "-r",
        "--inreads",
        metavar="READS",
        nargs="+",
        default=[],
        help="1 or more fastq files to extract reads from",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=False,
        metavar="OUTPUTSPREFIX",
        help="Output file prefix for reads. If disabled, then output it as all prefixes based on R1 or R2 or not",
    )
    return parser.parse_args(argv)
def open_file(file):
    """Open file, with support for gzip if the file ends with .gz."""
    if file.endswith(".gz"):
        return gzip.open(file, "rt")  # Open as text mode for .gz files
    else:
        return open(file, "r")

def parse_kraken2(k2_file, taxid_list):
    """Parse Kraken2 output and return a set of read IDs that match the taxid list."""
    read_ids = set()
    with open(k2_file, 'r') as infile:
        for line in infile:
            columns = line.strip().split()
            taxid_info = columns[-1]
            try:
                read_id = columns[1]  # Assuming second column is the read ID
                taxid_info_split = taxid_info.split(' ')
                taxid_list = set([t.split(":")[0] for t in taxid_info_split])
                if any(taxid in taxid_info for taxid in taxid_list):
                    read_ids.add(read_id)
            except Exception as e:
                print(f"Error parsing line: {e}")
    return read_ids
import re
# Remove using regex everything after the last occurrence of the extension for all removeexts from fastq_file
def remove_extension(fastq_file, removeexts):
    for ext in removeexts:
        fastq_file = re.sub(f"{ext}$", "", fastq_file)  # Use regex to match the extension at the end of the string
    return fastq_file
def match_files_base(inhandle, outhandle, read_ids):
    minreadlength = min(len(read_id) for read_id in read_ids)
    while True:
        header = inhandle.readline().strip()
        if not header:
            break  # EOF
        sequence = inhandle.readline().strip()
        plus = inhandle.readline().strip()
        quality = inhandle.readline().strip()

        # Extract the read ID (remove '@')
        raw_read_id = header[1:]
        raw_read_id = raw_read_id.split(' ')[0]  # Remove any metadata after space
        # Start with the full read ID, progressively remove characters until minreadlength
        match_found = False
        for length in range(len(raw_read_id), minreadlength - 1, -1):  # Start full, reduce to minreadlength
            trimmed_read_id = raw_read_id[:length]
            if trimmed_read_id in read_ids:
                # print(f"Matched read ID: {trimmed_read_id}")
                # Write the full read to the output file (header, sequence, +, quality)
                outhandle.write(f"{header}\n")
                outhandle.write(f"{sequence}\n")
                outhandle.write(f"{plus}\n")
                outhandle.write(f"{quality}\n")
                match_found = True
                break  # Exit the loop once a match is found

def match_files(in_handle, out_handle, read_ids):
    # iterate through the input FASTQ file with seqio.parse line by line
    # and filter out reads that are in the read_ids set

    # Write matching reads to output file as they are processed
    for record in SeqIO.parse(in_handle, "fastq"):
        # Extract the core part of the read ID before any space or metadata
        core_read_id = record.id.split()[0]
        core_read_id = core_read_id.split("/")[0]  # Remove any metadata after colon
        if core_read_id in read_ids:
            SeqIO.write(record, out_handle, "fastq")  # Write only matching reads


def match_reads(fastq_files, outprefix, read_ids):
    """Extract paired reads from fastq.gz files and write to new files dynamically."""
    # get minread length of read_ids
    # Handle paired-end reads
    if len(fastq_files) == 2:
        fastq1, fastq2 = fastq_files[0], fastq_files[1]

        # #### Prepare output files for R1 and R2
        out1 = f"{outprefix}.R1.removed.fastq"
        out2 = f"{outprefix}.R2.removed.fastq"
        with gzip.open(fastq1, "rt") as infq1, open(out1, "wt") as outgz1:
            match_files(infq1, outgz1, read_ids)
            # match_files_base(infq1, outgz1, read_ids)
        with gzip.open(fastq2, "rt") as infq2, open(out2, "wt") as outgz2:
            match_files(infq2, outgz2, read_ids)
            # match_files_base(infq2, outgz2, read_ids)
    else:
        fastq1 = fastq_files[0]
        out1 = f"{outprefix}.removed.fastq"
        with gzip.open(fastq1, "rt") as in_handle, open(out1, "wt") as out_handle:
            # Parse the input FASTQ and filter reads
            match_files(in_handle, out_handle, read_ids)

def main(argv=None):
    args = parse_args(argv)

    # Read the taxid list
    taxidfile = args.taxids
    taxid_list = []
    with open(taxidfile, 'r') as f:
        taxid_list = [line.strip() for line in f]

    # Parse Kraken2 classified reads file and get matching read IDs
    read_ids = parse_kraken2(args.k2, taxid_list)

    match_reads(args.inreads, args.prefix, read_ids)



if __name__ == "__main__":
    sys.exit(main())
