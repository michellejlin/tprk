import subprocess
import argparse
import sys
import os 
import regex
from itertools import chain
import numpy as np
from Bio.Seq import Seq
import pandas as pd

V1_list = list()	
V2_list = list()	
V3_list = list()	
V4_list = list()	
V5_list = list()	
V6_list = list()	
V7_list = list()	
V1_dna = list()
V2_dna = list()
V3_dna = list()
V4_dna = list()
V5_dna = list()
V6_dna = list()
V7_dna = list()

# Takes in the input file and goes through line by line, defining the read name and read sequence.
# Puts the read sequence through region_seq to match each read to the variable region/
def find_region(sample, input_format, is_pacbio, current_dir, strain_name):
	read_name = ""
	read_seq = ""
	region = ""
	sample_name = sample.split("." + input_format)[0]
	num_input_reads = 0

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
			# Illumina reads are read as one read = one count. PacBio reads are assumed to be RAD, and counts are
			# found after the underscore in the read name.
			if is_pacbio==False:
				# read_name = read_name.split('-')[0]
				read_count = 1
			# Includes the count in PacBio names after FAD/RAD, may be [1] or [2] depending on file
			else:
				read_count = int(float(read_name.split('_')[1].rstrip()))
			num_input_reads += read_count
		# Finds read sequence
		elif (line_num % line_block == 1):
			read_seq = str((Seq(line).reverse_complement()))
			# Matches each read seq to the region
			region_seq(read_seq, read_name, sample_name, is_pacbio, current_dir, num_input_reads)
		
	summary_table(num_input_reads, strain_name)

def summary_table(num_input_reads, strain_name):
	table = open(current_dir + "/" + strain_name + "_summary_statistics.csv", "w+")
	table.write("Sample,Total Input Reads,V1_Reads,V2_Reads,V3_Reads,V4_Reads,V5_Reads,V6_Reads,V7_Reads\n")
	to_write = strain_name + "," + str(num_input_reads) + "," + str(total_count(V1_list)) + "," + str(total_count(V2_list)) + "," + str(total_count(V3_list)) + "," + str(total_count(V4_list)) + "," + str(total_count(V5_list)) + "," + str(total_count(V6_list)) + "," + str(total_count(V7_list))
	table.write(to_write)

# Matches a read to a specified string of nucleotides, "primer", with less than 3 errors.
def fuzzy_match(read_seq, primer):
	# Finds exact match first. 
	exact_match = regex.search(primer,read_seq)
	if exact_match:
		return exact_match[0].rstrip()
	# If can't find exact match, searches for best match with less than 3 substitutions.
	else:
		fuzzy_match = regex.search(r"(?b)("+primer + "){s<=3}", read_seq)
		if fuzzy_match:
			return fuzzy_match[0]


# String matches a read to a variable region (V1-V7).
def region_seq(read_seq, read_name, sample_name, is_pacbio, current_dir, num_input_reads):
	# Different beginning portions for the variable regions.
	variable_regions = {
		"V1": "ATCAGTAGTAGTCTTAAATCC",
		#"V2": "CCGAACAAAATATCTCC",
		"V2": "AATATCTCCCCCCAATCCATA",
		# "V3": "GCCCCGACATCCCATAAGAT", #148B-148B2
		#"V3": "AGCACACAGACCCCAAAGCTT",
		"V3": "TCATACTCACCTTAGCCCCGAC",
		"V4": "AACAACGCATCTGCGCC",
		#"V5": "CCTTGGTTTCGAGCTT",
		"V5": "TCGAGCTTAATATAGGCAGC",
		"V6": "TTCCATACACCGGGAA",
		"V7": "CGGACTGACCACTACCCCACACTC"
	}

	# Loops through each variable region and its primer for each read.
	for region, primer in variable_regions.items():
		# If a read matches the primer with fuzzy matching less than 3 errors, then continue
		# Also reverses the read and tries that too
		read_seq_rev = str((Seq(read_seq).reverse_complement()))
		match_seq = fuzzy_match(read_seq, primer)
		match_seq_rev = fuzzy_match(read_seq_rev, primer)

		# Illumina reads are read as one read = one count. PacBio reads are assumed to be RAD, and counts are
		# found after the underscore in the read name.
		if is_pacbio==False:
			# read_name = read_name.split('-')[0]
			read_count = 1
		# Includes the count in PacBio names after FAD/RAD, may be [1] or [2] depending on file
		else:
			read_count = int(float(read_name.split('_')[1].rstrip()))

		if match_seq or match_seq_rev:
			end = 0
			if match_seq_rev:
				read_seq = read_seq_rev
				match_seq = match_seq_rev
			if (region == "V1"):
				start = (str.index(read_seq, match_seq) + 21)
				if fuzzy_match(read_seq, "CCAGGCCAGCTCCGCATA"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "CCAGGCCAGCTCCGCATA")))
					v_list = V1_list
					v_dna = V1_dna
			elif (region == "V2"):
				start = (str.index(read_seq, match_seq) + 21)
				if fuzzy_match(read_seq, "GTCGGTGTTAGACGCAAA"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "GTCGGTGTTAGACGCAAA")))
					v_list = V2_list
					v_dna = V2_dna
			elif (region == "V3"):
				start = (str.index(read_seq, match_seq) + 22)
				if fuzzy_match(read_seq, "GGAGTTGCCGGTGAGCTC"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "GGAGTTGCCGGTGAGCTC")))
					v_list = V3_list
					v_dna = V3_dna
			elif (region == "V4"):
				start = (str.index(read_seq, match_seq) + 15)
				if fuzzy_match(read_seq, "AGCAGCCAGAGCACACA"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "AGCAGCCAGAGCACACA")))
					v_list = V4_list
					v_dna = V4_dna
			elif (region == "V5"):
				# Find the beginning of the region following the primer
				start = (str.index(read_seq, match_seq) + 20)
				# Find the end of the region
				if fuzzy_match(read_seq, "CGATGCGAAATATCCTCC"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "CGATGCGAAATATCCTCC")))
					v_list = V5_list
					v_dna = V5_dna
			elif (region == "V6"):
				start = (str.index(read_seq, match_seq) + 16)
				if fuzzy_match(read_seq, "CATGTACGTACGCACATC"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "CATGTACGTACGCACATC")))
					v_list = V6_list
					v_dna = V6_dna
			elif (region == "V7"):
				start = (str.index(read_seq, match_seq) + 24)
				if fuzzy_match(read_seq, "CAAGTTTGCATACACTTT"):
					end = (str.index(read_seq, fuzzy_match(read_seq, "CAAGTTTGCATACACTTT")))	
					v_list = V7_list
					v_dna = V7_dna
			if end != 0:
				# Grabs what we now found as the region
				region_seq = read_seq[start:end]
				# Translates the sequence
				region_nuc = str((Seq(region_seq).reverse_complement()))
				region_seq = translate_nucs(region_nuc)
				all_assignments.write(sample_name + "," + read_name.rstrip()  + "," + region_seq + "," + 
					region_nuc + "," + region + "\n")

				# If the sequence is already in a global list of sequences for this region,
				# which is in [read sequence, count] format, increase the count by 1
				if (region_seq in chain(*v_list)):
					for sublist in v_list:
						if (sublist[0] == region_seq):
							sublist[1] = sublist[1] + read_count
				# Otherwise, add this sequence to the global list of sequences for this region
				else:
					v_list.append([region_seq, read_count])

				# Do the same thing for DNA.
				if (region_nuc in chain(*v_dna)):
					for sublist in v_dna:
						if (sublist[0] == region_nuc):
							sublist[1] = sublist[1] + read_count
				else:
					v_dna.append([region_nuc, read_count])

# Finds the total number of reads that mapped to a region.
def total_count(region_list):
	total = 0
	for read_seq, count in region_list:
		total = total + count
	return total

# Makes the final data table in csv format, with columns Region, Read, RelativeFreq, Count.
def make_table(strain_name, current_dir):
	Vlist_of_aas = [V1_list, V2_list, V3_list, V4_list, V5_list, V6_list, V7_list]
	Vlist_of_dna = [V1_dna, V2_dna, V3_dna, V4_dna, V5_dna, V6_dna, V7_dna]
	variable_regions = ["V1", "V2", "V3", "V4", "V5", "V6", "V7"]

	table = open(current_dir + "/" + strain_name + "_final_data.csv", "w+")
	table2 = open(current_dir + "/" + strain_name + "_over5count_final_dna_data.csv", "w+")
	table.write("Region,Read,RelativeFreq,Count" + "\n")
	table2.write("Region,Read,RelativeFreq,Count" + "\n")
	for index, v_list in enumerate(Vlist_of_aas):
		total = total_count(v_list)
		for read_seq, count in v_list:
			table.write(variable_regions[index] + "," + read_seq + "," + 
				str(((count / total) * 100)) + "," + str(count) + "\n")
	for index, v_list in enumerate(Vlist_of_dna):
		total = total_count(v_list)
		for read_seq, count in v_list:
			if(count>=5):
				table2.write(variable_regions[index] + "," + read_seq + "," + 
					str(((count / total) * 100)) + "," + str(count) + "\n")

	# Filters out the lines with greater than 1% to a separate _final_data_fitered.csv.
	#subprocess.call("awk -F\"[,|\\(]\" \'($3+0)>=1{print}\' " + strain_name + "_final_data.csv > " + strain_name + "_final_data_filtered.csv", shell=True)

# Translates a string of nucleotides into amino acids.
def translate_nucs(read_seq):
	coding_dna = Seq(read_seq)
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
	parser.add_argument('-s', '--strain_name', required=True,
		help='Specify what strain name.')
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
	file_name = args.strain_name

	# Right now the flag for input format is only to distinguish where the read name and read is.
	# Generally, fasta files are in a chunk of 2, fastq in a chunk of 4.
	if input_format in ["fasta", "fastq"]:
		print("Input format is " + input_format + ".")
	else:
		print("ERROR: Please specify either fasta or fastq for -i.")
		sys.exit()

	# Specifying whether the file is PacBio and Illumina changes how the program counts for each read.
	# Currently we assume PacBio sequences will be RAD-ified and in a specific format (i.e. >seq1_274) where
	# the reads have been clustered and the read count will be after one underscore.
	if is_pacbio and is_illumina:
		print("ERROR: -pacbio and -illumina are not compatible. Please select only one. ")
		sys.exit()
	if (is_pacbio == False) and (is_illumina==False):
		print("ERROR: Please select either -pacbio or -illumina.")
		sys.exit()

	# The all assignments file is a list of all the reads and what regions they mapped to, along with both
	# the nucleotide sequence and the amino acid sequence.
	all_assignments = open("all_assignments.csv", "w+")

	input_extension = input_format
	if input_format == "fasta" and is_pacbio:
		input_extension = "RAD.nolines.fix.fasta"

	# Extracts strain name from file name.
	strain_name = file_name.split("." + input_format)[0]
	# Get rid of stupidly long name in PacBio.
	if "RAD" in strain_name:
		strain_name = strain_name.split(".noprimers.filtered.RAD.nolines.fix.fasta")[0]
	
	# Checks if file exists already, and skips.
	file_tobemade = strain_name + "_final_data.csv"
	if os.path.isfile(file_tobemade):
		print(file_tobemade," already exists. Skipping making frequency tables...")
	else:
		print(file_tobemade,"")

		# Matches each read to a region and starts building a list.
		find_region(file_name, input_format, is_pacbio, current_dir, strain_name)
		
		# Builds the frequency final_table.csv for each file.
		make_table(strain_name, current_dir)
		V1_list = list()
		V2_list = list()	
		V3_list = list()	
		V4_list = list()	
		V5_list = list()	
		V6_list = list()	
		V7_list = list()
		V1_dna = list()
		V2_dna = list()
		V3_dna = list()
		V4_dna = list()
		V5_dna = list()
		V6_dna = list()
		V7_dna = list()