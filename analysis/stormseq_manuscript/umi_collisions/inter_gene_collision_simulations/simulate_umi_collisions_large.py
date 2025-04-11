#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import collections
import logging
import random
from scipy.stats import nbinom, poisson
import os

# Define the version number
VERSION = "0.0.2"

def parse_arguments():
    parser = argparse.ArgumentParser(description="UMI Distribution Simulation")
    parser.add_argument("--length_of_umi", type=int, default=10, help="Length of the UMI")
    parser.add_argument("--num_genes", type=int, default=100, help="Number of genes")
    parser.add_argument("--params_file", type=str, help="CSV file with empirical mu and theta values")
    parser.add_argument("--mu", type=float, help="Mean mu for negative binomial distribution")
    parser.add_argument("--theta", type=float, help="Dispersion theta for negative binomial distribution")
    parser.add_argument("--use_poisson", action='store_true', help="Use Poisson distribution instead of negative binomial if mu and theta are approximately equal")
    parser.add_argument("--gene_umi_counts_file", type=str, default="gene_umi_counts.txt", help="Output file for gene UMI counts")
    parser.add_argument("--umi_collisions_file", type=str, default="umi_inter_gene_collisions.txt", help="Output file for UMI inter-gene collisions")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    parser.add_argument("--log_file", type=str, default="simulation.log", help="Log file to redirect output")
    parser.add_argument("--umi_file", type=str, help="File containing UMIs to sample from")
    parser.add_argument("--version", action='version', version=f"%(prog)s {VERSION}", help="Show program's version number and exit")
    return parser.parse_args()

def setup_logging(log_file):
    logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def generate_all_umis(length):
    bases = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in itertools.product(bases, repeat=length)]

def load_umis_mmap(umi_file):
    # Count the number of lines in the file
    with open(umi_file, 'r') as f:
        num_lines = sum(1 for _ in f)
    # Create a memory-mapped array
    return np.memmap(umi_file, dtype='S16', mode='r', shape=(num_lines,))

def main():
    args = parse_arguments()
    setup_logging(args.log_file)
    
    length_of_umi = args.length_of_umi
    num_genes = args.num_genes
    params_file = args.params_file
    mu = args.mu
    theta = args.theta
    use_poisson = args.use_poisson
    gene_umi_counts_file = args.gene_umi_counts_file
    umi_collisions_file = args.umi_collisions_file
    seed = args.seed
    log_file = args.log_file
    umi_file = args.umi_file

    # Log the version number
    logging.info(f"Simulation Script Version: {VERSION}")

    # Log the parameters
    logging.info("Simulation Parameters:")
    logging.info(f"Length of UMI: {length_of_umi}")
    logging.info(f"Number of genes: {num_genes}")
    logging.info(f"Parameters file: {params_file}")
    logging.info(f"Mean (mu): {mu}")
    logging.info(f"Dispersion (theta): {theta}")
    logging.info(f"Use Poisson distribution: {use_poisson}")
    logging.info(f"Gene UMI counts file: {gene_umi_counts_file}")
    logging.info(f"UMI collisions file: {umi_collisions_file}")
    logging.info(f"Seed: {seed}")
    logging.info(f"Log file: {log_file}")
    logging.info(f"UMI file: {umi_file}")

    # If a seed is specified, use it. Otherwise, generate a random seed.
    if seed is not None:
        logging.info(f"Setting random seed to {seed}")
    else:
        seed = random.randint(0, 2**32 - 1)
        logging.info(f"No seed specified. Generated random seed: {seed}")
    np.random.seed(seed)
    random.seed(seed)

    if umi_file:
        logging.info(f"Using UMIs from file: {umi_file}")
        all_umis = load_umis_mmap(umi_file)
        total_umis = len(all_umis)
        logging.info(f"Total UMIs loaded: {total_umis}")
    else:
        logging.info("Generating all possible UMIs...")
        all_umis = generate_all_umis(length_of_umi)
        total_umis = len(all_umis)
        logging.info(f"Total possible UMIs generated: {total_umis}")

    logging.info("Loading parameters...")
    # Check if we are using empirical values or single mu and theta
    if params_file:
        # Load empirical mu and theta values
        params = pd.read_csv(params_file, index_col=0)
        mu_values = params['mean_mu'].values[:num_genes]
        theta_values = params['theta'].values[:num_genes]
    elif mu is not None and theta is not None:
        # Use provided single mu and theta for all genes
        mu_values = np.full(num_genes, mu)
        theta_values = np.full(num_genes, theta)
    else:
        raise ValueError("Either params_file or both mu and theta must be specified")

    logging.info("Sampling read depths for each gene...")
    # Sample read depths for each gene using empirical or provided mu and theta
    read_depths = []
    for mu, theta in zip(mu_values, theta_values):
        if use_poisson and np.isclose(mu, theta, rtol=0.1):
            # Use Poisson distribution
            read_depths.append(poisson.rvs(mu, size=1)[0])
        else:
            # Use Negative Binomial distribution
            p = theta / (mu + theta)
            read_depths.append(nbinom.rvs(mu, p, size=1)[0])
    logging.info("Read depths sampling complete.")

    logging.info("Allocating UMIs based on sampled read depths...")
    # Allocate UMIs based on the sampled read depths (sampling with replacement)
    allocated_umis = collections.defaultdict(list)
    for gene in range(num_genes):
        if gene % 100 == 0:
            logging.info(f"Allocating UMIs for gene {gene}...")
        count = read_depths[gene]
        if count > 0:
            sampled_umis = np.random.choice(all_umis, size=count, replace=True).astype(str)
            allocated_umis[gene].extend(sampled_umis)
    logging.info("UMI allocation complete.")

    logging.info("Counting overall UMI occurrences...")
    # Count overall UMI occurrences
    all_allocated_umis = [umi for umis in allocated_umis.values() for umi in umis]
    total_umis_allocated = len(all_allocated_umis)
    umi_counts = collections.Counter(all_allocated_umis)

    logging.info("Counting inter-gene collisions...")
    # Count inter-gene collisions
    umi_to_genes = collections.defaultdict(set)
    for gene, umis in allocated_umis.items():
        for umi in umis:
            umi_to_genes[umi].add(gene)

    inter_gene_collisions = {umi: len(genes) for umi, genes in umi_to_genes.items() if len(genes) > 1}

    # Calculate the inter-gene collision rate
    num_collisions = sum(len(genes) - 1 for genes in umi_to_genes.values() if len(genes) > 1)
    inter_gene_collision_rate = num_collisions / total_umis_allocated

    logging.info(f"Total UMIs allocated: {total_umis_allocated}")
    logging.info(f"Number of unique UMIs assigned to genes: {len(umi_counts)}")
    logging.info(f"Inter-gene collisions: {num_collisions}")
    logging.info(f"Inter-gene collision rate: {inter_gene_collision_rate:.5f}")

    logging.info("Writing gene UMI counts to file...")
    # Output gene UMI counts to a file
    with open(gene_umi_counts_file, 'w') as f:
        for gene, umis in allocated_umis.items():
            f.write(f"Gene {gene}: {len(umis)} UMIs\n")

    logging.info("Writing UMI inter-gene collisions to file...")
    # Output UMI inter-gene collisions to a file
    with open(umi_collisions_file, 'w') as f:
        for umi, collision_count in inter_gene_collisions.items():
            f.write(f"{umi}\t{collision_count}\n")

    logging.info("Simulation complete.")

if __name__ == "__main__":
    main()
