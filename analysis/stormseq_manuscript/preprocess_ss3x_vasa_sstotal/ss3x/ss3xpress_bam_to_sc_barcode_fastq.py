#!/usr/bin/env python3

import gzip
import argparse
from pathlib import Path

# Function to read barcodes from the external file
def read_barcodes(barcode_file_path):
    with open(barcode_file_path, 'r') as file:
        barcodes = [line.strip() for line in file]
    return barcodes

# Class to represent a FASTQ record
class FastqRecord:
    def __init__(self, id, sequence, quality):
        self.id = id
        self.sequence = sequence
        self.quality = quality

# Function to read a single FASTQ record from a gzipped file
def read_fastq_record(file):
    id = file.readline().decode().strip()
    if not id:
        return None
    sequence = file.readline().decode().strip()
    file.readline()  # Skip the '+' line
    quality = file.readline().decode().strip()
    return FastqRecord(id, sequence, quality)

# Function to process a batch of barcodes
def process_batch(barcodes_batch, i1_path, i2_path, r1_path, r2_path, output_path):
    file_streams = {}

    with gzip.open(i1_path, 'rb') as i1_file, \
         gzip.open(i2_path, 'rb') as i2_file, \
         gzip.open(r1_path, 'rb') as r1_file, \
         gzip.open(r2_path, 'rb') as r2_file:

        while True:
            i1_record = read_fastq_record(i1_file)
            i2_record = read_fastq_record(i2_file)
            r1_record = read_fastq_record(r1_file)
            r2_record = read_fastq_record(r2_file)

            if not i1_record or not i2_record or not r1_record or not r2_record:
                break

            for barcode in barcodes_batch:
                if barcode in i1_record.sequence:
                    if barcode not in file_streams:
                        file_streams[barcode] = {
                            'I1': open(output_path / f"{barcode}_I1.fastq", 'w'),
                            'I2': open(output_path / f"{barcode}_I2.fastq", 'w'),
                            'R1': open(output_path / f"{barcode}_R1.fastq", 'w'),
                            'R2': open(output_path / f"{barcode}_R2.fastq", 'w')
                        }

                    streams = file_streams[barcode]
                    streams['I1'].write(f"{i1_record.id}\n{i1_record.sequence}\n+\n{i1_record.quality}\n")
                    streams['I2'].write(f"{i2_record.id}\n{i2_record.sequence}\n+\n{i2_record.quality}\n")
                    streams['R1'].write(f"{r1_record.id}\n{r1_record.sequence}\n+\n{r1_record.quality}\n")
                    streams['R2'].write(f"{r2_record.id}\n{r2_record.sequence}\n+\n{r2_record.quality}\n")
                    break

    # Close all file streams
    for barcode_streams in file_streams.values():
        for stream in barcode_streams.values():
            stream.close()

# Parse command line arguments
parser = argparse.ArgumentParser(description='Process FASTQ files based on barcodes.')
parser.add_argument('-b', '--barcode_file', required=True, help='Path to the barcode file')
parser.add_argument('-i1', '--i1_path', required=True, help='Path to I1 fastq.gz file')
parser.add_argument('-i2', '--i2_path', required=True, help='Path to I2 fastq.gz file')
parser.add_argument('-r1', '--r1_path', required=True, help='Path to R1 fastq.gz file')
parser.add_argument('-r2', '--r2_path', required=True, help='Path to R2 fastq.gz file')
parser.add_argument('-o', '--output_path', required=True, help='Output directory path')
parser.add_argument('-n', '--batch_size', type=int, required=True, help='Number of barcodes per batch')
args = parser.parse_args()

barcodes = read_barcodes(args.barcode_file)
output_path = Path(args.output_path)
output_path.mkdir(parents=True, exist_ok=True)

# Process barcodes in batches
for i in range(0, len(barcodes), args.batch_size):
    batch = barcodes[i:i + args.batch_size]
    process_batch(batch, args.i1_path, args.i2_path, args.r1_path, args.r2_path, output_path)
