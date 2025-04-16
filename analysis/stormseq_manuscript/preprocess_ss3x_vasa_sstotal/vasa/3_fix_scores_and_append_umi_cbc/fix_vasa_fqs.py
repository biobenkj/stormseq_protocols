#!/usr/bin/env python3

import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import argparse

def modify_fastq_record(description):
    tags = {tag.split(':')[0]: tag.split(':')[1] for tag in description.split(';')[1:]}
    tags['QT'] = tags['QT'].upper().replace('O', '/')
    tags['RQ'] = tags['RQ'].upper().replace('O', '/')
    new_tags = ';'.join([f"{key}:{tags[key]}" for key in sorted(tags.keys())])
    return description.split(';')[0] + ';' + new_tags, tags

def modify_sequence_and_quality(record, tags):
    # Generate the new sequence
    modified_seq_str = tags['RX'] + tags['CB'] + str(record.seq)
    modified_seq = Seq(modified_seq_str)
    
    # Convert quality scores to Phred scores, prepend RQ and QT qualities
    qt_rq_quality_scores = [ord(c)-33 for c in tags['RQ'] + tags['QT']]
    modified_quality_scores = qt_rq_quality_scores + record.letter_annotations['phred_quality']
    
    # Create a new SeqRecord with the modified sequence and reassign quality scores
    modified_record = SeqRecord(modified_seq,
                                id=record.id,
                                name=record.name,
                                description=record.description,
                                letter_annotations={'phred_quality': modified_quality_scores})
    return modified_record

def process_fastq(fastq_gz_path, output_gz_path):
    with gzip.open(fastq_gz_path, 'rt') as input_gz, gzip.open(output_gz_path, 'wt') as output_gz:
        for record in SeqIO.parse(input_gz, "fastq"):
            modified_description, tags = modify_fastq_record(record.description)
            record.description = modified_description
            record.id = modified_description.split()[0]
            
            record = modify_sequence_and_quality(record, tags)
            
            SeqIO.write(record, output_gz, "fastq")

def main():
    parser = argparse.ArgumentParser(description="Modify FASTQ records.")
    parser.add_argument("input_file", help="Input gzipped FASTQ file")
    parser.add_argument("output_file", help="Output gzipped FASTQ file")
    args = parser.parse_args()

    process_fastq(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

