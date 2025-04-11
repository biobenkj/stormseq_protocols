#!/usr/bin/env python3

import pysam
from collections import defaultdict
import argparse

def calculate_umi_collisions(bam_file_path, output_file_path=None, filter_dash=False, collision_bam_path=None):
    # Open the BAM file
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")

    # Create a new BAM file for writing collisions if specified
    collision_bam_file = None
    if collision_bam_path:
        collision_bam_file = pysam.AlignmentFile(collision_bam_path, "wb", header=bam_file.header)

    # Dictionary to store UMI and gene associations
    umi_gene_dict = defaultdict(set)
    umi_read_dict = defaultdict(set)  # Use set to store unique (gene, chromosome, start, end, read_name) tuples

    # Iterate over each read in the BAM file
    for read in bam_file:
        # Skip unmapped reads, secondary alignments, and supplementary alignments
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Extract UMI (UB tag, or UR tag if UB is '-')
        try:
            umi = read.get_tag("UB")
            if umi == "-":
                umi = read.get_tag("UR")
            gene = read.get_tag("GX")
            start = read.reference_start
            end = read.reference_end
            read_name = read.query_name

            # Add the gene and coordinates to the set associated with this UMI
            umi_gene_dict[umi].add(gene)
            umi_read_dict[umi].add((gene, read.reference_name, start, end, read_name))
        except KeyError:
            # Skip reads without UB, UR, or GX tags
            continue

    bam_file.close()

    # Identify inter-gene UMI collisions
    inter_gene_collisions = {umi: genes for umi, genes in umi_gene_dict.items() if len(genes) > 1}

    # Filter out UMIs associated with "-" entries if requested
    if filter_dash:
        inter_gene_collisions = {umi: genes for umi, genes in inter_gene_collisions.items() if '-' not in genes}

    # Print the total number of inter-gene UMI collisions
    print("Number of inter-gene UMI collisions:", len(inter_gene_collisions))

    # Write the results to the output file if specified
    if output_file_path:
        with open(output_file_path, 'w') as outfile:
            outfile.write("UMI\tGene\tChromosome\tStart\tEnd\tReadName\n")
            for umi in inter_gene_collisions:
                for gene, chrom, start, end, read_name in umi_read_dict[umi]:
                    outfile.write(f"{umi}\t{gene}\t{chrom}\t{start}\t{end}\t{read_name}\n")

    # Reopen the BAM file to fetch and write collision reads if specified
    if collision_bam_file:
        bam_file = pysam.AlignmentFile(bam_file_path, "rb")
        for umi in inter_gene_collisions.keys():
            for gene, chrom, start, end, read_name in umi_read_dict[umi]:
                for read in bam_file.fetch(chrom, start, end):
                    try:
                        read_umi = read.get_tag("UB")
                        if read_umi == "-":
                            read_umi = read.get_tag("UR")
                        if read_umi == umi:
                            collision_bam_file.write(read)
                    except KeyError:
                        continue
        collision_bam_file.close()

def main():
    parser = argparse.ArgumentParser(description='Calculate inter-gene UMI collisions from STARsolo BAM files.')
    parser.add_argument('bam_file', type=str, help='Path to the STARsolo BAM file')
    parser.add_argument('-o', '--output', type=str, help='Path to the output text file', default=None)
    parser.add_argument('--filter-dash', action='store_true', help='Filter out UMIs associated with "-" entries')
    parser.add_argument('--collision-bam', type=str, help='Path to the output BAM file for collision reads', default=None)

    args = parser.parse_args()
    calculate_umi_collisions(args.bam_file, args.output, args.filter_dash, args.collision_bam)

if __name__ == '__main__':
    main()
