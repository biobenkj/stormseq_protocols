#!/usr/bin/env python3

import gzip
from Bio import SeqIO
import argparse

def filter_reads(input_fastq, output_fastq):
    with gzip.open(input_fastq, 'rt') as infile, gzip.open(output_fastq, 'wt') as outfile:
        for record in SeqIO.parse(infile, "fastq"):
            if record.seq.endswith("TTT"):
                SeqIO.write(record, outfile, "fastq")

def main():
    parser = argparse.ArgumentParser(description="Filter reads from a gzipped FASTQ file that end with 'TTT'.")
    parser.add_argument("input_fastq", type=str, help="Input gzipped FASTQ file")
    parser.add_argument("output_fastq", type=str, help="Output gzipped FASTQ file")

    args = parser.parse_args()

    filter_reads(args.input_fastq, args.output_fastq)

if __name__ == "__main__":
    main()
