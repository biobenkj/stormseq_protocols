#!/usr/bin/env python3
import gzip
import argparse

def process_ss3_reads(read1_file, read2_file, output_prefix):
    pattern = 'ATTGCGCAATG'
    with gzip.open(read1_file, 'rt') as r1, gzip.open(read2_file, 'rt') as r2, \
         gzip.open(f'{output_prefix}_five_prime_r1.fastq.gz', 'wt') as r1_five, \
         gzip.open(f'{output_prefix}_five_prime_r2.fastq.gz', 'wt') as r2_five, \
         gzip.open(f'{output_prefix}_internal_reads_r1.fastq.gz', 'wt') as r1_internal, \
         gzip.open(f'{output_prefix}_internal_reads_r2.fastq.gz', 'wt') as r2_internal:

        while True:
            # Read four lines at a time (FASTQ format)
            read1 = [r1.readline().strip() for _ in range(4)]
            read2 = [r2.readline().strip() for _ in range(4)]

            # Break if end of file is reached
            if not read1[0] or not read2[0]:
                break

            # Check if read1 starts with the pattern
            if read1[1].startswith(pattern):
                for line in read1:
                    r1_five.write(line + '\n')
                for line in read2:
                    r2_five.write(line + '\n')
            else:
                for line in read1:
                    r1_internal.write(line + '\n')
                for line in read2:
                    r2_internal.write(line + '\n')

def main():
    parser = argparse.ArgumentParser(description='Process paired-end FASTQ files.')
    parser.add_argument('read1_file', type=str, help='Path to read 1 FASTQ file (gzipped)')
    parser.add_argument('read2_file', type=str, help='Path to read 2 FASTQ file (gzipped)')
    parser.add_argument('output_prefix', type=str, help='Prefix for the output files')

    args = parser.parse_args()

    process_ss3_reads(args.read1_file, args.read2_file, args.output_prefix)

if __name__ == '__main__':
    main()
