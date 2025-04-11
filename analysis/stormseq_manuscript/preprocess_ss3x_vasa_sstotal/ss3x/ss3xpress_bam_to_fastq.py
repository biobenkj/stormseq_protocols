#!/usr/bin/env python3

import pysam
import argparse
import gzip

def convert_ss3_ubam_to_fastq(bam_file, fastq_prefix, ub_tag_length):
    # Set the check_sq to False as this is an unmapped bam
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        with gzip.open(f"{fastq_prefix}_I1.fastq.gz", "wt") as i1_file, \
             gzip.open(f"{fastq_prefix}_I2.fastq.gz", "wt") as i2_file, \
             gzip.open(f"{fastq_prefix}_R1.fastq.gz", "wt") as r1_file, \
             gzip.open(f"{fastq_prefix}_R2.fastq.gz", "wt") as r2_file:

            for read in bam.fetch(until_eof=True):
                # Extract tags
                bc = read.get_tag("BC")
                qb = read.get_tag("QB")
                ub = read.get_tag("UB")
                qu = read.get_tag("QU")

                # Determine read name
                read_name = f"@{read.query_name}"

                # Process R1 and R2 FASTQ files
                if read.is_read1:
                    seq, qual = read.query_sequence, read.qual
                    # Check UB tag length and prepend sequence if necessary
                    if len(ub) == ub_tag_length:
                        seq = "ATTGCGCAATG" + ub + "GGG" + seq
                        qual = "IIIIIIIIIII" + qu + "III" + qual
                    r1_file.write(f"{read_name}\n{seq}\n+\n{qual}\n")
                    # Write barcode info only for first read in pair
                    i1_file.write(f"{read_name}\n{bc}\n+\n{qb}\n")
                elif read.is_read2:
                    r2_file.write(f"{read_name}\n{read.query_sequence}\n+\n{read.qual}\n")
                    # Only write to I2 for read2 to avoid duplication
                    i2_file.write(f"{read_name}\n{bc}\n+\n{qb}\n")

def main():
    parser = argparse.ArgumentParser(description='Convert BAM to FASTQ')
    parser.add_argument('bam_file', type=str, help='Input BAM file')
    parser.add_argument('fastq_prefix', type=str, help='Prefix for the output FASTQ files')
    parser.add_argument('ub_tag_length', type=int, help='Expected length of the UB BAM tag')

    args = parser.parse_args()
    convert_ss3_ubam_to_fastq(args.bam_file, args.fastq_prefix, args.ub_tag_length)

if __name__ == "__main__":
    main()
