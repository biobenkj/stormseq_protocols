# TEProf3: TE-derived Promoter Finder 3
# Author: Holden Liang
# Contributor: Ben K Johnson
# -----------------------------------------------------------------------------
# Functions related to processing bam files and performing de novo assembly

import argparse
import time
import subprocess
import os
os.environ['NUMEXPR_MAX_THREADS'] = '4'
os.environ['NUMEXPR_NUM_THREADS'] = '2'
import numexpr as ne
import multiprocess as mp
from itertools import repeat

import misc

def move_gtf(input_dataset):
	"""
	Move GTF files to assembled folder if the user 
	provided GTF files directly (assemble_mode=0).
	"""
	for sample in input_dataset:
		try:
			old_path = input_dataset[sample]['gtf']
			new_path = "assembled/" + input_dataset[sample]['gtf']
			subprocess.run(f"mv {old_path} {new_path}", shell=True, executable='/bin/bash')
			input_dataset[sample]['gtf'] = new_path
		except:
			misc.print_time("No GTF file is provided. (1) Set --assemblemode to run de novo assembly, or (2) Check if you put GTFs in your sample_manifest.txt.")
			exit()
	return input_dataset

def start_assemble(input_dataset, flags):
	"""
	Run different assembly commands based on the user's choice:
	  0 = no assembly (use provided GTF)
	  1 = short-read only
	  2 = long-read only
	  3 = short+long in a single hybrid run (--mix)
	  4 = short & long read assembled separately in two passes
	"""
	assemble_mode = flags.assemblemode 
	misc.print_time("Run de novo transcript assembly (mode="+assemble_mode+")")

	if assemble_mode == "0":
		misc.print_time("You chose NOT to run de novo assembly, using user-provided GTF. Moving GTFs to assembled/ folder.")
		input_dataset = move_gtf(input_dataset)

	elif assemble_mode == "1":
		misc.print_time("You chose to run StringTie on your SHORT-READ data only.")
		input_dataset = stringtie_short_read(input_dataset, flags)

	elif assemble_mode == "2":
		misc.print_time("You chose to run StringTie on your LONG-READ data only.")
		input_dataset = stringtie_long_read(input_dataset, flags)

	elif assemble_mode == "3":
		misc.print_time("You chose to run a HYBRID assembly (short + long together in one step).")
		input_dataset = stringtie_hybrid_read(input_dataset, flags)

	elif assemble_mode == "4":
		misc.print_time("You chose to run SHORT and LONG reads in separate passes.")
		# We'll do short-read assembly for samples that have short reads,
		# and long-read assembly for samples that have long reads.
		input_dataset = stringtie_separate_short_and_long(input_dataset, flags)

	else:
		misc.print_time(f"Unrecognized assemblemode: {assemble_mode}")
		exit()

	misc.print_time("Assembly step completed.")
	return input_dataset


def stringtie_short_read(input_dataset, flags):
	"""
	Run StringTie on short-read data for each sample. 
	Use per-sample strand code if provided. 
	"""
	assemble_thread = flags.assemblethread
	assemble_length = flags.assemblelength
	assemble_junctionread = flags.assemblejunctionread
	assemble_samplenumber = flags.assemblesamplenumber

	commands = []
	with open('./assembled/assemble_commands.txt', 'w') as f:
		for sample in input_dataset:
			if "bam_SR" not in input_dataset[sample]:
				continue  # skip samples that lack short reads

			sample_strand_code = input_dataset[sample].get("SR_strand", flags.assemblestrand)
			if sample_strand_code == '0':
				assemble_strand = ''
			elif sample_strand_code == '1':
				assemble_strand = '--rf'
			elif sample_strand_code == '2':
				assemble_strand = '--fr'
			else:
				misc.print_time(f"Unrecognized strand code {sample_strand_code} for sample {sample}")
				exit()

			input_bam = input_dataset[sample]["bam_SR"]
			out_gtf   = "assembled/" + os.path.basename(input_bam).replace("bam","gtf")

			command = ' '.join([
				"samtools view -q 255 -h", 
				input_bam, 
				"| stringtie",
				"-",
				"-o", out_gtf,
				"-j", assemble_junctionread,
				"-p", assemble_thread,
				"-m", assemble_length,
				assemble_strand
			])

			commands.append(command)
			f.write(command + "\n")

			input_dataset[sample]["gtf"] = out_gtf

	# Run commands in parallel
	if commands:
		with mp.Pool(assemble_samplenumber) as pool:
			pool.map(misc.run_bash_command, commands)
	else:
		misc.print_time("No short-read samples found to assemble.")
	return input_dataset

def stringtie_long_read(input_dataset, flags):
	"""
	Run StringTie on long-read data for each sample. 
	Use per-sample strand code if provided.
	"""
	assemble_thread = flags.assemblethread
	assemble_length = flags.assemblelength
	assemble_junctionread = flags.assemblejunctionread
	assemble_samplenumber = flags.assemblesamplenumber

	commands = []
	with open('./assembled/assemble_commands.txt', 'w') as f:
		for sample in input_dataset:
			if "bam_LR" not in input_dataset[sample]:
				continue  # skip samples that lack long reads

			sample_strand_code = input_dataset[sample].get("LR_strand", flags.assemblestrand)
			if sample_strand_code == '0':
				assemble_strand = ''
			elif sample_strand_code == '1':
				assemble_strand = '--rf'
			elif sample_strand_code == '2':
				assemble_strand = '--fr'
			else:
				misc.print_time(f"Unrecognized strand code {sample_strand_code} for sample {sample}")
				exit()

			input_LR_bam = input_dataset[sample]["bam_LR"]
			out_gtf = "assembled/" + os.path.basename(input_LR_bam).replace("bam","gtf")

			command = ' '.join([
				"stringtie -L",
				"-o", out_gtf,
				"-j", assemble_junctionread,
				"-p", assemble_thread,
				"-m", assemble_length,
				assemble_strand,
				input_LR_bam
			])

			commands.append(command)
			f.write(command + "\n")

			input_dataset[sample]["gtf"] = out_gtf

	if commands:
		with mp.Pool(assemble_samplenumber) as pool:
			pool.map(misc.run_bash_command, commands)
	else:
		misc.print_time("No long-read samples found to assemble.")
	return input_dataset

def stringtie_hybrid_read(input_dataset, flags):
	"""
	Run StringTie --mix for short+long reads in one command per sample.
	If the short-read strand code != long-read strand code, 
	force unstranded and warn the user.
	"""
	assemble_thread = flags.assemblethread
	assemble_length = flags.assemblelength
	assemble_junctionread = flags.assemblejunctionread
	assemble_samplenumber = flags.assemblesamplenumber

	def code_to_stringtie_arg(code):
		if code == '0':
			return ''
		elif code == '1':
			return '--rf'
		elif code == '2':
			return '--fr'
		else:
			misc.print_time(f"Unrecognized strand code: {code}")
			exit()

	commands = []
	with open('./assembled/assemble_commands.txt', 'w') as f:
		for sample in input_dataset:
			# Need both short and long for a hybrid run
			if "bam_SR" not in input_dataset[sample] or "bam_LR" not in input_dataset[sample]:
				continue

			sr_code = input_dataset[sample].get("SR_strand", flags.assemblestrand)
			lr_code = input_dataset[sample].get("LR_strand", flags.assemblestrand)

			if sr_code != lr_code:
				misc.print_time(
					f"Warning: Sample {sample} short-read strand ({sr_code}) differs "
					f"from long-read strand ({lr_code}). Forcing unstranded assembly."
				)
				assemble_strand_arg = ''
			else:
				assemble_strand_arg = code_to_stringtie_arg(sr_code)

			input_SR_bam = input_dataset[sample]["bam_SR"]
			input_LR_bam = input_dataset[sample]["bam_LR"]

			out_gtf = "assembled/" + os.path.basename(input_SR_bam).replace("bam","gtf")
			cmd_list = [
				"stringtie",
				"--mix",
				"-o", out_gtf,
				"-j", assemble_junctionread,
				"-p", assemble_thread,
				"-m", assemble_length,
				assemble_strand_arg,
				input_SR_bam,
				input_LR_bam
			]
			command = ' '.join(cmd_list)
			f.write(command + "\n")
			commands.append(command)

			input_dataset[sample]["gtf"] = out_gtf

	if commands:
		with mp.Pool(assemble_samplenumber) as pool:
			pool.map(misc.run_bash_command, commands)
	else:
		misc.print_time("No sample had both short and long reads for a hybrid assembly.")

	return input_dataset

def stringtie_separate_short_and_long(input_dataset, flags):
	"""
	New mode (4): assemble short reads and long reads 
	in separate passes for all samples that have them.
	Essentially do 'stringtie_short_read' then 'stringtie_long_read'.
	"""
	misc.print_time("First, assemble short reads (if any).")
	input_dataset = stringtie_short_read(input_dataset, flags)
	misc.print_time("Now, assemble long reads (if any).")
	input_dataset = stringtie_long_read(input_dataset, flags)
	return input_dataset


def check_strand_info(input_dataset, flags):
	"""
	Check if the GTF files produced (or provided) 
	have actual strand annotation (+/-) 
	"""
	misc.print_time("Check if there's strand information in the GTF file(s)")
	input_gtf_files = []
	for sample in input_dataset:
		input_gtf_files.append(input_dataset[sample]['gtf'])

	no_strand_flags = []
	for gtf in input_gtf_files:
		with open(gtf, "r") as input_file:
			no_strand_flag = 1
			for line in input_file:
				if line.startswith("#"):
					continue
				entry = line.strip().split("\t")
				# Skip alt scaffolds or mitochondria if you do so
				if "GL" not in entry[0] and "KI" not in entry[0] and "chrM" not in entry[0]:
					if entry[2] == "transcript":
						strand = entry[6]
						if strand in ["+", "-"]:
							no_strand_flag = 0
							break
			if no_strand_flag == 1:
				misc.print_time("This GTF file doesn't have strand info: " + gtf)
			no_strand_flags.append(no_strand_flag)
	return max(no_strand_flags)
