# Syphilis Project
import subprocess
import argparse
import sys
import os 
import regex
from itertools import chain
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd

V1_list = list()	
V2_list = list()	
V3_list = list()	
V4_list = list()	
V5_list = list()	
V6_list = list()	
V7_list = list()	
V1_count = 0	
V2_count = 0	
V3_count = 0	
V4_count = 0	
V5_count = 0	
V6_count = 0	
V7_count = 0

# Takes in the input file and goes through line by line, defining the read name and read sequence.
# Puts the read sequence through region_seq to match each read to the variable region/
def find_region(sample, input_format, is_pacbio, current_dir):
	read_name = ""
	read_seq = ""
	region = ""
	sample_name = sample.split("." + input_format)[0]

	# Fastq files are in blocks of 4: with first line being name, second line being sequence,
	# whereas fasta files are in blocks of 2.
	if input_format == "fasta":
		line_block = 2
	else:
		line_block = 4

	for line_num, line in enumerate(open(current_dir + "/" + sample)):
		# Finds read name
		if (line_num % line_block == 0):
			read_name = line
		# Finds read sequence
		elif (line_num % line_block == 1):
			read_seq = str((Seq(line, generic_dna).reverse_complement()))
			# Matches each read seq to the region
			region_seq(read_seq, read_name, sample_name, is_pacbio, current_dir)

# Matches a read to a specified string of nucleotides, "primer", with less than 3 errors.
def fuzzy_match(read_seq, primer):
	match = regex.findall(primer + "{e<=3}", read_seq)
	if match:
		return match[0]
	else:
		return False

# String matches a read to a variable region (V1-V7).
def region_seq(read_seq, read_name, sample_name, is_pacbio, current_dir):
	# Different beginning portions for the variable regions.
	variable_regions = {
		"V1": "ATCAGTAGTAGTCTTAAATCC",
		"V2": "CCGAACAAAATATCTCC",
		"V3": "GCCCCGACATCCCATAAGAT",
		"V4": "AACAACGCATCTGCGCC",
		"V5": "CCTTGGTTTCGAGCTT",
		"V6": "CATACACCGGGAAGTT",
		"V7": "TACCCCACACTC"
	}

	# Loops through each variable region and its primer for each read.
	for region, primer in variable_regions.items():
		# If a read matches the primer with fuzzy matching less than 3 errors, then continue
		# Also reverses the read and tries that too

		read_seq_rev = str((Seq(read_seq, generic_dna).reverse_complement()))
		match_seq = regex.findall(primer + "{e<=3}", read_seq)
		match_seq_rev = regex.findall(primer + "{e<=3}", read_seq_rev)

		if is_pacbio==False:
			# read_name = read_name.split('-')[0]
			read_count = 1
		# Includes the count in PacBio names after FAD/RAD, may be [1] or [2] depending on file
		else:
			read_count = int(float(read_name.split('_')[1].rstrip()))

		if match_seq or match_seq_rev:
			if match_seq_rev:
				read_seq = read_seq_rev
				match_seq = match_seq_rev
			if (region == "V1"):
				start = (str.index(read_seq, match_seq[0]) + 21)
				if fuzzy_match(read_seq, "CCAGGCCAGCTCCGCA"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "CCAGGCCAGCTCCGCA")))
					# Grabs what we now found as the region
					region_seq = read_seq[start:end]
					# Translates the sequence
					region_nuc = str((Seq(region_seq, generic_dna).reverse_complement()))
					region_seq = translate_nucs(region_nuc)
					all_assignments.write(read_name.rstrip()  + "," + region_seq + "," + region_nuc + ",V1"+ "\n")

					# If the sequence is already in a global list of sequences for this region,
					# which is in [read sequence, count] format, increase the count by 1
					if (region_seq in chain(*V1_list)):
						for sublist in V1_list:
							if (sublist[0] == region_seq):
								sublist[1] = sublist[1] + read_count
					# Otherwise, add this sequence to the global list of sequences for this region
					else:
						V1_list.append([region_seq, read_count])

			elif (region == "V2"):
				start = (str.index(read_seq, match_seq[0]) + 17)
				if fuzzy_match(read_seq, "GTCGGTGTTAGACGCAA"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "GTCGGTGTTAGACGCAA")))
					region_seq = read_seq[start:end]
					region_nuc = str((Seq(region_seq, generic_dna).reverse_complement()))
					region_seq = translate_nucs(region_nuc)
					all_assignments.write(read_name.rstrip()  + "," + region_seq + "," + region_nuc + ",V2"+ "\n")
					if (region_seq in chain(*V2_list)):
						for sublist in V2_list:
							if (sublist[0] == region_seq):
								sublist[1] = sublist[1] + read_count
					else:
						V2_list.append([region_seq, read_count])

			elif (region == "V3"):
				start = (str.index(read_seq, match_seq[0]) + 20)
				if fuzzy_match(read_seq, "GGAGTTGCCGGT"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "GGAGTTGCCGGT")))
					region_seq = read_seq[start:end]
					region_nuc = str((Seq(region_seq, generic_dna).reverse_complement()))
					region_seq = translate_nucs(region_nuc)
					all_assignments.write(read_name.rstrip()  + "," + region_seq + "," + region_nuc + ",V3"+ "\n")
					if (region_seq in chain(*V3_list)):
						for sublist in V3_list:
							if (sublist[0] == region_seq):
								sublist[1] = sublist[1] + read_count
					else:
						V3_list.append([region_seq, read_count])

			elif (region == "V4"):
				start = (str.index(read_seq, match_seq[0]) + 15)
				if fuzzy_match(read_seq, "AGCAGCCAGAGCACAC"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "AGCAGCCAGAGCACAC")))
					region_seq = read_seq[start:end]
					region_nuc = str((Seq(region_seq, generic_dna).reverse_complement()))
					region_seq = translate_nucs(region_nuc)
					all_assignments.write(read_name.rstrip()  + "," + region_seq + "," + region_nuc + ",V4"+ "\n")
					if (region_seq in chain(*V4_list)):
						for sublist in V4_list:
							if (sublist[0] == region_seq):
								sublist[1] = sublist[1] + read_count
					else:
						V4_list.append([region_seq, read_count])

			elif (region == "V5"):
				# Find the beginning of the region following the primer
				start = (str.index(read_seq, match_seq[0]) + 16)
				# Find the end of the region
				if fuzzy_match(read_seq, "CGATGCGAAATA"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "CGATGCGAAATA")))
					region_seq = read_seq[start:end]
					region_nuc = str((Seq(region_seq, generic_dna).reverse_complement()))
					region_seq = translate_nucs(region_nuc)
					all_assignments.write(read_name.rstrip()  + "," + region_seq + "," + region_nuc + ",V5"+ "\n")
					if (region_seq in chain(*V5_list)):
						for sublist in V5_list:
							if (sublist[0] == region_seq):
								sublist[1] = sublist[1] + read_count
					else:
						V5_list.append([region_seq, read_count])

			elif (region == "V6"):
				start = (str.index(read_seq, match_seq[0]) + 16)
				if fuzzy_match(read_seq, "CATGTACGTACG"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "CATGTACGTACG")))
					region_seq = read_seq[start:end]
					region_nuc = str((Seq(region_seq, generic_dna).reverse_complement()))
					region_seq = translate_nucs(region_nuc)
					all_assignments.write(read_name.rstrip()  + "," + region_seq + "," + region_nuc + ",V6"+ "\n")
					if (region_seq in chain(*V6_list)):
						for sublist in V6_list:
							if (sublist[0] == region_seq):
								sublist[1] = sublist[1] + read_count
					else:
						V6_list.append([region_seq, read_count])

			elif (region == "V7"):
				start = (str.index(read_seq, match_seq[0]) + 12)
				if fuzzy_match(read_seq, "CAAGTTTGCATACACTT"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "CAAGTTTGCATACACTT")))	
					region_seq = read_seq[start:end]
					region_nuc = str((Seq(region_seq, generic_dna).reverse_complement()))
					region_seq = translate_nucs(region_nuc)
					all_assignments.write(read_name.rstrip()  + "," + region_seq + "," + region_nuc + ",V7"+ "\n")
					if (region_seq in chain(*V7_list)):
						for sublist in V7_list:
							if (sublist[0] == region_seq):
								sublist[1] = sublist[1] + read_count
					else:
						V7_list.append([region_seq, read_count])

# Finds the total number of reads that mapped to a region.
def total_count(region_list):
	total = 0
	for read_seq, count in region_list:
		total = total + count
	return total

# Makes the final data table in csv format, with columns Region, Read, RelativeFreq, Count.
def make_table(strain_name, current_dir):
	table = open(current_dir + "/" + strain_name + "_final_data.csv", "w+")
	table.write("Region,Read,RelativeFreq, Count" + "\n")
	for read_seq, count in V1_list:
		total = total_count(V1_list)
		table.write("V1," + read_seq + "," + str(((count / total) * 100)) + "," + str(count) + "\n")
	for read_seq, count in V2_list:
		total = total_count(V2_list)
		table.write("V2," + read_seq + "," + str(((count / total) * 100)) + "," + str(count) + "\n")
	for read_seq, count in V3_list:
		total = total_count(V3_list)
		table.write("V3," + read_seq + "," + str(((count / total) * 100)) + "," + str(count) + "\n")

	for read_seq, count in V4_list:
		total = total_count(V4_list)
		table.write("V4," + read_seq + "," + str(((count / total) * 100)) + "," + str(count) + "\n")

	for read_seq, count in V5_list:
		total = total_count(V5_list)
		table.write("V5," + read_seq + "," + str(((count / total) * 100)) + "," + str(count) + "\n")

	for read_seq, count in V6_list:
		total = total_count(V6_list)
		table.write("V6," + read_seq + "," + str(((count / total) * 100)) + "," + str(count) + "\n")

	for read_seq, count in V7_list:
		total = total_count(V7_list)
		table.write("V7," + read_seq + "," + str(((count / total) * 100)) + "," + str(count) + "\n")

	# Filters out the lines with greater than 1% to a separate _final_data_fitered.csv.
	#subprocess.call("awk -F\"[,|\\(]\" \'($3+0)>=1{print}\' " + strain_name + "_final_data.csv > " + strain_name + "_final_data_filtered.csv", shell=True)

# Translates a string of nucleotides into amino acids.
def translate_nucs(read_seq):
	coding_dna = Seq(read_seq, generic_dna)
	translation = str(coding_dna.translate())
	return translation

if __name__ == '__main__': 
	parser = argparse.ArgumentParser(description='tprK project. Example usage: syph.py -i fasta -p pacbio')
	parser.add_argument('-i', '--input_format', required=True,
		help='Specify either fasta or fastq. Will take all the files in this folder with the specified extension. ')
	parser.add_argument('-pacbio', action='store_true', help='Write this flag to specify this file is a pacbio file. '
		'Do not use with -illumina') 
	parser.add_argument('-illumina', action='store_true', help='Write this flag to specify this file is an illumina file. ' 
		'Do not use with -pacbio')
	parser.add_argument('-d', '--directory', help='Pass directory (used when passed in from within R.')
	
	# Checks for argument sanity.
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)
	
	current_dir = args.directory

	input_format = args.input_format.lower()
	is_pacbio = args.pacbio
	is_illumina = args.illumina

	if input_format in ["fasta", "fastq"]:
		print("Input format is " + input_format + ". Reading all " + input_format + " files in this folder. ")
	else:
		print("ERROR: Please specify either fasta or fastq for -i.")
		sys.exit()

	if is_pacbio and is_illumina:
		print("ERROR: -pacbio and -illumina are not compatible. Please select only one. ")
		sys.exit()
	if (is_pacbio == False) and (is_illumina==False):
		print("ERROR: Please select either -pacbio or -illumina.")
		sys.exit()

	all_assignments = open("all_assignments.csv", "w+")

	for file in os.listdir(current_dir):
		if file.endswith(input_format):
			find_region(file, input_format, is_pacbio, current_dir)

			strain_name = file.split("." + input_format)[0]
			make_table(strain_name, current_dir)
			V1_list = list()
			V2_list = list()	
			V3_list = list()	
			V4_list = list()	
			V5_list = list()	
			V6_list = list()	
			V7_list = list()	
			V1_count = 0	
			V2_count = 0	
			V3_count = 0	
			V4_count = 0	
			V5_count = 0	
			V6_count = 0	
			V7_count = 0
	
	# Call the visualizer.
	# subprocess.call('python3 ' + dir_path + '/syph_visualizer.py ' + current_dir + '/final_data.csv', shell=True)
	# subprocess.call('python3 ' + dir_path + '/syph_network.py ' + current_dir + '/all_assignments.csv ' + current_dir + '/final_data.csv', shell=True)

# samtools view 148B-TPRK3.bam|awk '{print$3 "\t" $4 "\t" $4+length($10)-1}' > pos.txt
