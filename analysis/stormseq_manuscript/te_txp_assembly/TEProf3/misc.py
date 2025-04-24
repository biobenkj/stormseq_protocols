# TEProf3: TE-derived Promoter Finder 3
# Author: Holden Liang
# Contributor: Ben K Johnson
# -----------------------------------------------------------------------------
# Functions related to general utility, argument parsing, and manifest parsing

import argparse
import time
import subprocess
import os
os.environ['NUMEXPR_MAX_THREADS'] = '4'
os.environ['NUMEXPR_NUM_THREADS'] = '2'
import numexpr as ne
import multiprocess as mp
from multiprocessing import Pool
from itertools import repeat
import sys
import pandas as pd
from tqdm import tqdm
import glob
import numpy as np
import ast
import json

def get_args():
    """ Fetches the arguments for the program """
    program_desc = """TEProf3 takes aligned files (bam files) and/or assembled data (gtf files)
    and provides you a list of TE-derived transcripts with expression across samples."""
    parser = argparse.ArgumentParser(description=program_desc)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    ## for debugging
    parser.add_argument("--test", dest='test',
                        help="testing mode", action='store_true')

    ## bonus functions
    parser.add_argument("--teprof2", dest='teprof2',
                        help="directly use teprof2 output for translation", type=str, default="no")
    parser.add_argument("--gdcjson", dest='gdcjson',
                        help="json file downloaded from GDC for splice junction tsv files to generate sample_manifest.txt file", type=str, default="no")
    parser.add_argument("--gtexjson", dest='gtexjson',
                        help="json file downloaded from GTEx for splice junction tsv files to generate sample_manifest.txt file", type=str, default="no")

    ## for identifying neoantigen candidate 
    parser.add_argument("--split", dest='split',
                        help="run blastp on sbatch cluster like htcf, put the fasta file here", type=str, default="no")
    parser.add_argument("--blastp", dest='blastp',
                        help="run blastp on sbatch cluster like htcf, put the fasta file here", type=str, default="no")
    parser.add_argument("--blastpshort", dest='blastpshort',
                        help="run blastp-short on sbatch cluster like htcf, put the fasta file here", type=str, default="no")
    parser.add_argument("--blastpdatabase", dest='blastpdatabase',
                        help="database that blastp will run against (e.g., nr or gene)", type=str, default="nr")
    parser.add_argument("--taxiddb", dest='taxiddb',
                        help="folder to taxid file", type=str, 
                        default=os.path.dirname(os.path.realpath(sys.argv[0])) + "/../reference/blast_database_nr/")
    parser.add_argument("--blastjobnum", dest='blastjobnum',
                        help="number of blast jobs run simultaneously", type=str, default="50")
    parser.add_argument("--classify", dest='classify',
                        help="path to teprof3_protein_information.tsv table from teprof3 --translation output", type=str, default="no")  
    parser.add_argument("--detectedprotein", dest='detectedprotein',
                        help="process BLAST and BLAT result to get a list of detected TE-derived proteins", action="store_true") 
    parser.add_argument("--parsenetmhcpan", dest='parsenetmhcpan',
                        help="parse output from netMHCpan", type=str, default="no") 
    parser.add_argument("--neoantigen", dest='neoantigen',
                        help="get list of neoantigen candidates", action="store_true") 

    ## prepare reference files
    parser.add_argument("--repeatmasker", dest='repeatmasker',
                        help="repeatmasker file, only provide it if you are making a new reference", type=str, default="no")
    parser.add_argument("--geneannotation", dest='geneannotation',
                        help="gene annotation file, only provide it if you are making a new reference. GENCODE gtf is recommended", type=str, default="no")
    parser.add_argument("--geneannotationfortranslation", dest='geneannotationfortranslation',
                        help="gene annotation file for translation, only provide it if you are making a new reference. GENCODE gtf is recommended", action="store_true")
    parser.add_argument("--hervannotation", dest='hervannotation',
                        help="herv annotation file", type=str, default="no")
    parser.add_argument("--geneprotein", dest='geneprotein',
                        help="generate protein sequence fasta file from gene annotation", action="store_true")

    ## general parameters
    parser.add_argument("-f", "--manifest", dest="manifest",
                        help="Dataset manifest file: sample name, short/long read, bam/gtf file (ended with bam or gtf). (tab-delimited). Optionally add a 4th column for strand code (0,1,2).",
                        type=str)
    parser.add_argument("-ki", "--keepintermediate", dest="keepintermediate",
                        help="Keep intermediate files", action="store_true")
    parser.add_argument("-s", "--samplenumber", dest="samplenumber",
                        help="number of samples to be processed together at the same time (default 10)", type=int, default=10)
    parser.add_argument("-g", "--guided", dest="guided",
                        help="run teprof3 in guided mode...", type=str, default="no")
    parser.add_argument("-rs", "--reset", dest="reset",
                        help="reset the environment from a previous run", action="store_true")
    parser.add_argument("-v", "--version", dest="version",
                        help="Print version of teprof3", action="store_true")

    ## assemble parameter
    parser.add_argument("-am", "--assemblemode", dest='assemblemode',
                        help=(
                            "how to run transcript de novo assembly:\n"
                            "0 = none (use provided GTF),\n"
                            "1 = short-read only,\n"
                            "2 = long-read only,\n"
                            "3 = short+long in a single hybrid run (--mix),\n"
                            "4 = short & long assembled separately.\n"
                            "(default 0)"
                        ),
                        type=str, default="0")
    parser.add_argument("-at", "--assemblethread", dest='assemblethread',
                        help="number of threads used for assembly (default 4)", type=str, default="4")
    parser.add_argument("-al", "--assemblelength", dest='assemblelength',
                        help="minimum transcript length for stringtie assembly (default 200)", type=str, default="200")
    parser.add_argument("-as", "--assemblesamplenumber", dest='assemblesamplenumber',
                        help="number of samples to be processed in parallel (default 10)", type=int, default=10)
    parser.add_argument("-ast", "--assemblestrand", dest='assemblestrand',
                        help="default strandedness (0=unstranded, 1=--rf, 2=--fr). (default 0)", type=str, default="0")
    parser.add_argument("-aj", "--assemblejunctionread", dest='assemblejunctionread',
                        help="minimum junction coverage (default 1)", type=str, default="1")

    ## process assemble parameter
    parser.add_argument("-pt", "--processtpm", dest='processtpm',
                        help="TPM cutoff to filter te-derived transcripts right after assembly (default 0.5)",
                        type=float, default=0.5)
    parser.add_argument("-ps", "--processsamplenumber", dest='processsamplenumber',
                        help="number of samples to process at the same time (default 10)",
                        type=int, default=10)
    parser.add_argument("-ptn", "--processtranscriptnumber", dest='processtranscriptnumber',
                        help="only include samples with > X TE-derived transcripts in sample. (default 100)",
                        type=int, default=100)

    ## transcript filtering parameter
    parser.add_argument("-fm", "--filtermode", dest='filtermode',
                        help="how to filter TE-derived transcripts (1=only short read, 2=short+long) (default 1)",
                        type=str, default="1")
    parser.add_argument("-fs", "--filtersamplenumber", dest='filtersamplenumber',
                        help="number of samples to process at once (default 10)", type=int, default=10)
    parser.add_argument("-ft", "--filterthread", dest='filterthread',
                        help="threads for filtering transcripts (default 10)", type=int, default=10)
    parser.add_argument("-fa", "--filterannotated", dest='filterannotated',
                        help="keep annotated transcripts. Default is to exclude them.",
                        action="store_false")
    parser.add_argument("-fi", "--filterintronretention", dest='filterintronretention',
                        help="cutoff for intron retention filtering. (default 3)",
                        type=int, default=3)
    parser.add_argument("-fmo", "--filtermonoexon", dest='filtermonoexon',
                        help="exclude mono-exonic transcripts (default false).",
                        action='store_true')
    parser.add_argument("-fmot", "--filtermonoexontpm", dest='filtermonoexontpm',
                        help="TPM cutoff for mono-exonic transcripts (default 1).",
                        type=float, default=1)
    parser.add_argument("-fdm", "--filterdownstreammate", dest='filterdownstreammate',
                        help="# reads capturing splicing junction to downstream exon (default 2).",
                        type=int, default=2)
    parser.add_argument("-fr", "--filterratio", dest='filterratio',
                        help="ratio used in chimeric mate filtering (default 0.5).",
                        type=float, default=0.5)
    parser.add_argument("-fncf", "--filternochimericfilter", dest='filternochimericfilter',
                        help="disable chimeric mate filter", action='store_true')
    parser.add_argument("-fljt", "--longreadsjtolerance", dest='longreadsjtolerance',
                        help="SJ tolerance for LR vs stringtie (default 3)",
                        type=int, default=3)
    parser.add_argument("-flmt", "--longreadmtolerance", dest='longreadmtolerance',
                        help="tolerance for mono-exonic check in LR data (default 50)",
                        type=int, default=50)

    ## mega assembly parameter
    parser.add_argument("-tt", "--tacothread", dest='tacothread',
                        help="number of threads for TACO (default 10)",
                        type=str, default="10")
    parser.add_argument("-tto", "--tacotolerance", dest='tacotolerance',
                        help="edge tolerance for TACO (default 3)", type=int, default=3)

    ## quantification parameter
    parser.add_argument("-qm", "--quanmode", dest='quanmode',
                        help="1=short read, 2=SJ only (default 1)",
                        type=str, default="1")
    parser.add_argument("-qs", "--quansamplenumber", dest='quansamplenumber',
                        help="samples in parallel for stringtie quant (default 10)",
                        type=int, default=10)
    parser.add_argument("-qsc", "--quansamplenumbercon", dest='quansamplenumbercon',
                        help="samples in parallel for concatenating stringtie quant output (default 50)",
                        type=int, default=50)
    parser.add_argument("-qnp", "--quannoprepde", dest='quannoprepde',
                        help="skip prepDE.py for raw counts (default False)",
                        action='store_true')
    parser.add_argument("-qng", "--quannogene", dest='quannogene',
                        help="collect TPM for canonical genes (default True)",
                        action='store_true')
    parser.add_argument("-ql", "--quanreadlength", dest='quanreadlength',
                        help="read length of libraries (default 75)",
                        type=str, default="75")

    ## in silico translation parameter
    parser.add_argument("-ti", "--translation", dest='translation',
                        help="provide GTF or table for translation", type=str, default="no")
    parser.add_argument("-tm", "--translationmode", dest='translationmode',
                        help="translation mode (1 or 2, default 1)", type=str, default="1")
    parser.add_argument("-tl", "--translationlength", dest='translationlength',
                        help="minimum length of translated protein (default 20)",
                        type=int, default=20)
    parser.add_argument("-tg", "--translationgenome", dest="translationgenome",
                        help="reference genome name for translation (default hg38)",
                        type=str, default="hg38")

    args = parser.parse_args()
    return args

def print_version(version):
    """ Save version of each tool used in the pipeline """
    subprocess.run("echo \"TEProf3 version:\" > debug/teprof3_version_info.txt", shell=True, executable='/bin/bash')
    subprocess.run("echo \""+version+"\" >> debug/teprof3_version_info.txt", shell=True, executable='/bin/bash')
    subprocess.run("echo \"Stringtie\" >> debug/teprof3_version_info.txt", shell=True, executable='/bin/bash')
    subprocess.run("stringtie --version >> debug/teprof3_version_info.txt", shell=True, executable='/bin/bash')


def parse_manifest(input_file, default_strand):
    """
    Parse the input manifest. 
    Allow either 3 or 4 columns:
     1) sample name
     2) library type (short/long/gtf/SJ)
     3) path to file
     4) OPTIONAL strand code (0, 1, or 2). 
        If missing, fallback to 'default_strand' from flags.assemblestrand.
    Stores:
     - short reads => 'bam_SR', 'SR_strand'
     - long  reads => 'bam_LR', 'LR_strand'
     - gtf => 'gtf'
     - SJ => 'SJ'
    Returns a dict: input_dict[sampleName][key] = value
    """

    if not os.path.isfile(input_file):
        print_time(f"Can't find {input_file}. Please check.")
        return 1

    input_dict = {}
    with open(input_file, 'r') as file:
        for line in file:
            entry = line.strip().split('\t')
            # Skip empty lines
            if entry == ['']:
                continue

            # Accept lines with 3 or 4 columns
            if len(entry) == 3 or len(entry) == 4:
                sampleName  = entry[0]
                libraryType = entry[1]
                pathToFile  = entry[2]

                # If 4th col is present, use it; else fallback
                if len(entry) == 4:
                    strandCode = entry[3]
                else:
                    strandCode = default_strand  # fallback

                # Initialize dict for sample if needed
                if sampleName not in input_dict:
                    input_dict[sampleName] = {}

                if libraryType == "long":
                    input_dict[sampleName]["bam_LR"]    = pathToFile
                    input_dict[sampleName]["LR_strand"] = strandCode
                elif libraryType == "short":
                    input_dict[sampleName]["bam_SR"]    = pathToFile
                    input_dict[sampleName]["SR_strand"] = strandCode
                elif libraryType == "gtf":
                    input_dict[sampleName]["gtf"]       = pathToFile
                elif libraryType == "SJ":
                    input_dict[sampleName]["SJ"]        = pathToFile
                else:
                    print_time(f"Unrecognized libraryType '{libraryType}' in manifest. Expect short, long, gtf, or SJ.")
                    return 1
            else:
                print_time("Please provide 3 or 4 columns per line in the manifest file.")
                return 1

    return input_dict

def print_dataset(input_dataset):
    """ print input dataset for confirmation and record """
    for sample in input_dataset:
        for key in input_dataset[sample]:
            print('\t'.join([sample, key, str(input_dataset[sample][key])]))

def move_file(input_dataset, flags):
    """
    Move bam files and GTF files (if present) to subfolders 
    for a cleaner folder space.
    """
    print_time("Move files to folders")
    for sample in input_dataset:

        # Move long-read BAM if present
        if "bam_LR" in input_dataset[sample]:
            lr_bam = input_dataset[sample]["bam_LR"]
            try:
                subprocess.run(f"mv {lr_bam}.bai bam/{lr_bam}.bai", shell=True, executable='/bin/bash')
            except:
                pass  # no index for LR possibly

            try:
                subprocess.run(f"mv {lr_bam} bam/{lr_bam}", shell=True, executable='/bin/bash')
                input_dataset[sample]["bam_LR"] = "bam/" + lr_bam
            except:
                print_time(f"No long read bam file found for {sample}")

        # Move short-read BAM if present
        if "bam_SR" in input_dataset[sample]:
            sr_bam = input_dataset[sample]["bam_SR"]
            try:
                subprocess.run(f"mv {sr_bam}.bai bam/{sr_bam}.bai", shell=True, executable='/bin/bash')
            except:
                if flags.quanmode == "1":
                    print_time(f"No index file for short read bam files is provided for {sample}")
                    return 1

            try:
                subprocess.run(f"mv {sr_bam} bam/{sr_bam}", shell=True, executable='/bin/bash')
                input_dataset[sample]["bam_SR"] = "bam/" + sr_bam
            except:
                if flags.quanmode == "1":
                    print_time(f"No short read bam file found for {sample}")
                    return 1

        # Move SJ if present
        if "SJ" in input_dataset[sample]:
            sj_file = input_dataset[sample]["SJ"]
            try:
                subprocess.run(f"mv {sj_file} bam/{sj_file}", shell=True, executable='/bin/bash')
                input_dataset[sample]["SJ"] = "bam/" + sj_file
            except:
                print_time(f"No SJ file is provided for {sample}")
                return 1

        # Move GTF if present
        if "gtf" in input_dataset[sample]:
            gtf_file = input_dataset[sample]["gtf"]
            try:
                # Move the sample's GTF to the assembled folder
                subprocess.run(f"mv {gtf_file} assembled/{gtf_file}", shell=True, executable='/bin/bash')
                input_dataset[sample]["gtf"] = "assembled/" + gtf_file
            except:
                print_time(f"Could not find GTF file {gtf_file} for {sample}")
                return 1

    return input_dataset   

def print_time(input_string):
    """ Print message with time stamp """
    program_start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    print("[ "+program_start_time+" ] "+input_string )

def run_bash_command(command):
    subprocess.run(command, shell=True, executable='/bin/bash')

def gtf_to_refbed(input_gtf_file):
    """ convert stringtie gtf file to refbed file for WashU genome browser visualization """
    output_refbed_file = input_gtf_file.replace("gtf", "refbed")
    input_data = {}
    with open(input_gtf_file, 'r') as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue
            entry = line.strip('\n').split('\t')
            detail_info = entry[8]
            transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
            try:
                gene_name = detail_info.split("gene_name \"")[1].split("\";")[0]
            except:
                gene_name = transcript_id
            transcript_type = "coding"
            if entry[2] == "transcript":
                input_data[transcript_id] = [[
                    entry[0], 
                    str(int(entry[3])-1), 
                    entry[4],
                    str(int(entry[3])-1),
                    entry[4],
                    entry[6],
                    gene_name,
                    transcript_id,
                    transcript_type
                ],[],[],[detail_info]]
            elif entry[2] == "exon":
                input_data[transcript_id][1].append(int(entry[3])-1)
                input_data[transcript_id][2].append(int(entry[4]))

    with open(output_refbed_file, 'w') as output_file:
        for tx_id in input_data:
            item = input_data[tx_id]
            exon_starts = sorted(item[1])
            exon_ends   = sorted(item[2])
            output_line = item[0] + [
                ','.join(str(x) for x in exon_starts),
                ','.join(str(x) for x in exon_ends)
            ] + item[3]
            output_file.write('\t'.join(output_line) + '\n')

    run_bash_command(f"sort -k1,1 -k2,2n {output_refbed_file} | bgzip > {output_refbed_file}.sorted.gz")
    run_bash_command(f"tabix -p bed {output_refbed_file}.sorted.gz")
    run_bash_command(f"rm {output_refbed_file}")

def all_gtf_to_refbed(input_dataset, flags):
    gtf_files = []
    if flags.guided == "no":
        for sample in input_dataset:
            gtf_files.append(input_dataset[sample]['gtf'])
        gtf_files += glob.glob("assembled/*.TE.gtf")
    else:
        gtf_files.append("assembled/TE_transcript_consensus.gtf")
    gtf_files.append("assembled/TE_transcript_consensus.filtered.gtf")

    with mp.Pool(20) as pool:
        pool.map(gtf_to_refbed, gtf_files)

def process_SJ_gdc_json_to_prep_sample_manifest(input_json_file):
    with open("sample_manifest.txt","w") as output_file:
        input_data = json.load(open(input_json_file,"r"))
        for sample in tqdm(input_data, desc="Extracting information:"):
            if sample["data_format"] == "TSV":
                case_id = sample["associated_entities"][0]["case_id"]
                entity_id = sample["associated_entities"][0]["entity_id"]
                bam_file = sample["analysis"]["input_files"][0]["file_name"]
                sj_file  = sample["file_name"].replace(".gz","")
                # Write short read line
                output_file.write('\t'.join([f"{case_id}_{entity_id}", "short", bam_file]) + "\n")
                # Write gtf line
                output_file.write('\t'.join([f"{case_id}_{entity_id}", "gtf",  bam_file.replace("bam","gtf")]) + "\n")
                # Write SJ line
                output_file.write('\t'.join([f"{case_id}_{entity_id}", "SJ",   sj_file]) + "\n")

def process_SJ_gtex_json_to_prep_sample_manifest(input_json_file):
    with open("sample_manifest.txt","w") as output_file:
        input_data = json.load(open(input_json_file,"r"))
        for sample in tqdm(input_data, desc="Extracting information:"):
            obj_id = sample["object_id"].split("/")[1]
            file_name = sample["file_name"]
            output_file.write('\t'.join([obj_id, "SJ", file_name])+"\n")

def reset_folder(flags):
    input_file = flags.manifest
    if not os.path.isfile(input_file):
        print_time("Can't find " + input_file + ". Please check.")
        exit()
    files_in_assembled = []
    files_in_bam = []
    with open(input_file, 'r') as file:
        for line in file:
            entry = line.strip().split('\t')
            if len(entry) >= 3:
                lib_type = entry[1]
                fpath    = entry[2]
                if lib_type in ["long", "short"]:
                    files_in_bam.append(fpath)
                    files_in_bam.append(fpath+".bai")
                elif lib_type == "gtf":
                    files_in_assembled.append(fpath)
                elif lib_type == "SJ":
                    files_in_bam.append(fpath)
            elif entry == ['']:
                continue
            else:
                print_time("please put at least 3 arguments in each row in the manifest file")
                exit()

    print_time("reset files from a failed run")
    for f in files_in_assembled:
        run_bash_command(f"mv ./assembled/{f} ./{f}")
    for f in files_in_bam:
        run_bash_command(f"mv ./bam/{f} ./{f}")

    run_bash_command("rm -rf bam assembled debug")
    return None

def cleanup_folder():
    run_bash_command("mkdir assembled/intermediate_files; mv assembled/transcript_count_matrix.csv assembled/intermediate_files/transcript_count_matrix.csv; mv assembled/taco_command.txt assembled/intermediate_files/taco_command.txt; mv assembled/sample_list_for_prepDE.txt assembled/intermediate_files/sample_list_for_prepDE.txt")
    run_bash_command("mv assembled/quantification_command.txt assembled/intermediate_files/quantification_command.txt; mv assembled/TE_transcript_consensus.gtf assembled/intermediate_files/TE_transcript_consensus.gtf; mv assembled/TACO_output assembled/intermediate_files/TACO_output")
    run_bash_command("mv assembled/gtf_to_merge.txt assembled/intermediate_files/gtf_to_merge.txt; mv assembled/gene_count_matrix.csv assembled/intermediate_files/gene_count_matrix.csv; mv assembled/*quantification.gtf assembled/intermediate_files/; mv assembled/*stats assembled/intermediate_files/; mv assembled/*TE.gtf assembled/intermediate_files/")
    run_bash_command("mv assembled/TE_transcript_consensus.filtered.gtf assembled/teprof3_output_TE_transcript_consensus.gtf; mv assembled/TE_transcript_consensus.filtered.refbed.sorted.gz assembled/teprof3_output_TE_transcript_consensus.refbed.sorted.gz; mv assembled/TE_transcript_consensus.filtered.refbed.sorted.gz.tbi assembled/teprof3_output_TE_transcript_consensus.refbed.sorted.gz.tbi")
    run_bash_command("gzip assembled/teprof3_output_quantification.tsv")
    run_bash_command("gzip assembled/teprof3_output_quantification.TE.tsv")


def summarize_stas(input_dataset):
    """
    Summarize .stats.txt files, but ONLY for short-read samples.
    This avoids reading stats from purely long-read samples
    which may have inconsistent stats lines or no stats at all.
    """

    gtf_to_sample_dict = {}
    input_gtf_files = []

    # Only gather GTFs for samples that have short-read data
    for sample in input_dataset:
        # Check if the sample has "bam_SR" (meaning short reads)
        if "bam_SR" in input_dataset[sample]:
            # This sample has short-read data => we expect a .stats.txt
            gtf_file = input_dataset[sample]["gtf"]
            gtf_to_sample_dict[gtf_file] = sample
            input_gtf_files.append(gtf_file)

    stats_dict = {}
    for input_gtf_file in input_gtf_files:
        sample_name = gtf_to_sample_dict[input_gtf_file]
        stats_dict[sample_name] = []
        header = []

        # Read the .stats.txt lines for this GTF
        stats_file_path = input_gtf_file + ".stats.txt"
        try:
            with open(stats_file_path, "r") as stats_file:
                for line in stats_file:
                    entry = line.strip("\n").split("\t")
                    # entry[1] => the step name/description
                    # entry[2] => the numeric value
                    stats_dict[sample_name].append(entry[2])
                    header.append(entry[1])
        except FileNotFoundError:
            # Maybe this sample has no .stats.txt, skip or log a warning
            print(f"Warning: no stats file found for {input_gtf_file}")
            continue

    # Build a DataFrame from stats_dict
    # stats_dict keys: sample names => list of values
    # 'header' was re-set for each file, so we only have the final one
    # If you expect the same lines across all short-read samples,
    # you can just take the header from the last iteration
    stats_pandas = pd.DataFrame.from_dict(stats_dict, orient='index')

    # If at least one sample produced lines, we can set columns
    if not stats_pandas.empty:
        stats_pandas.columns = header
        stats_pandas = stats_pandas.transpose()
        stats_pandas.to_csv("assembled/teprof3_output_transcript_statistic.tsv", sep="\t")

        # Remove the .stats.txt files if desired
        run_bash_command("rm assembled/*.gtf.stats.txt")
    else:
        print("No short-read stats were found or no consistent stats to summarize.")

        
