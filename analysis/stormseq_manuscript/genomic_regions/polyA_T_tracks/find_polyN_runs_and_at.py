#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from multiprocessing import Pool, cpu_count

def parse_arguments():
    parser = argparse.ArgumentParser(description='Identify polyN tracks using a greedy algorithm with consecutive mismatches.')
    parser.add_argument('-i', '--input', required=True, help='Input genome FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    parser.add_argument('-n', '--nucleotide', required=True, choices=['A', 'T', 'C', 'G'], help='Nucleotide to search for')
    parser.add_argument('-l', '--min_bases', type=int, default=5, help='Minimum number of starting nucleotides required in a run (default: 5)')
    parser.add_argument('-c', '--max_consecutive_mismatches', type=int, default=1, help='Maximum number of consecutive mismatches allowed before ending the run (default: 1)')
    parser.add_argument('--filter_at_content', action='store_true', help='Filter runs to keep only those with at least 60% AT content')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads (default: 1)')
    return parser.parse_args()

def calculate_at_content(seq):
    """Calculate the AT content of a given sequence."""
    at_bases = seq.count('A') + seq.count('T')
    total_bases = len(seq)
    if total_bases == 0:
        return 0
    return at_bases / total_bases

def process_chromosome(args):
    record, nucleotide, min_bases, max_consecutive_mismatches, filter_at_content = args
    seq = str(record.seq).upper()
    chrom = record.id
    seq_length = len(seq)
    results = []

    in_run = False
    start = 0
    last_starting_base_position = 0
    num_starting_bases = 0
    consecutive_mismatches = 0

    for i in range(seq_length):
        base = seq[i]
        if base == nucleotide:
            if not in_run:
                # Start a new run
                in_run = True
                start = i
                num_starting_bases = 1
                consecutive_mismatches = 0
            else:
                # Continue the run
                num_starting_bases += 1
                consecutive_mismatches = 0
            last_starting_base_position = i  # Update last starting base position
        else:
            if in_run:
                consecutive_mismatches += 1
                if consecutive_mismatches > max_consecutive_mismatches:
                    # End the run at the last starting nucleotide
                    end = last_starting_base_position + 1
                    if num_starting_bases >= min_bases:
                        run_seq = seq[start:end]
                        if not filter_at_content or calculate_at_content(run_seq) >= 0.6:
                            results.append({
                                'chrom': chrom,
                                'start': start,
                                'end': end,
                                'nucleotide': nucleotide,
                                'num_starting_bases': num_starting_bases
                            })
                    in_run = False
                    num_starting_bases = 0
                    consecutive_mismatches = 0
                else:
                    # Continue the run despite the mismatch
                    continue
            else:
                # Not in a run, do nothing
                continue

    # Handle run at the end of the sequence
    if in_run:
        # Adjust end position based on mismatches at the end
        if consecutive_mismatches >= 1:
            end = last_starting_base_position + 1  # End at last starting nucleotide
        else:
            end = seq_length
        if num_starting_bases >= min_bases:
            run_seq = seq[start:end]
            if not filter_at_content or calculate_at_content(run_seq) >= 0.6:
                results.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'nucleotide': nucleotide,
                    'num_starting_bases': num_starting_bases
                })

    return results

def main():
    args = parse_arguments()

    input_file = args.input
    output_file = args.output
    nucleotide = args.nucleotide.upper()
    min_bases = args.min_bases
    max_consecutive_mismatches = args.max_consecutive_mismatches
    filter_at_content = args.filter_at_content
    num_threads = args.threads

    # Read sequences
    records = list(SeqIO.parse(input_file, "fasta"))

    # Prepare arguments for multiprocessing
    pool_args = [(record, nucleotide, min_bases, max_consecutive_mismatches, filter_at_content) for record in records]

    # Process chromosomes
    if num_threads > 1:
        with Pool(processes=num_threads) as pool:
            all_results = pool.map(process_chromosome, pool_args)
    else:
        all_results = [process_chromosome(arg) for arg in pool_args]

    # Flatten results
    flattened_results = [run for result in all_results for run in result]

    # Write to output file
    with open(output_file, 'w') as f:
        for run in flattened_results:
            chrom = run['chrom']
            start = run['start']
            end = run['end']
            nucleotide = run['nucleotide']
            num_starting_bases = run['num_starting_bases']
            f.write(f"{chrom}\t{start}\t{end}\t{nucleotide}_run_{num_starting_bases}\n")

if __name__ == "__main__":
    main()
