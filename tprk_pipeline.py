# tprK pipeline
import subprocess
import argparse
import sys
import os 
import regex

# Filters allreads.csv and outputs allreads_filtered.csv based on flags -f and -c.
# (Default relative freq is 0.2 and count is 5.) Thrown away reads are put into allreads_thrownout.csv.
def filter_table(relative_freq_filter, count_filter, table):
	line_num = 0
	table_csv = open(table, "r")
	table_lines = table_csv.readlines()
	table_reads_thrownout = open(table.split(".csv")[0] + "_thrownout.csv", "w+")

	print ("Filtering " + table + " by count > " + str(count_filter) + " and relative frequency > " + str(relative_freq_filter) + "...")
	for line in table_lines:
		line = line.rstrip()
		# Ignores the first line and initializes files.
		if line_num == 0:
			first_line = line
			table_filtered = open(table.split(".csv")[0] + "_filtered.csv", "w+")
			table_filtered.write(line + "\n")
			table_reads_thrownout.write(line + "\n")
		# Filters lines based on specified values. Checks that each sample passes both the relative frequency filter
		# and the count filter. Writes passed lines to alldata_filtered.csv, and thrown away lines to alldata_thrownout.csv.
		else:
			line_parts = line.split(",")
			relative_freq_check = False
			count_check = False
			for index, part in enumerate(line_parts):
				# Even indexes should be relative frequencies
				if ((index % 2) == 0) and index != 0 and index != 1 and part != "NA" and float(part) >= relative_freq_filter:
					relative_freq_check = True
				# Odd indexes should be counts
				if ((index % 2 == 1)) and index != 0 and index != 1 and part != "NA" and float(part) >= count_filter:
					count_check = True
				# For allreads_filtered.csv, we want to change anything under the count filter to NA.
				if ("allreads" in table) and ((index % 2 == 1)) and index != 0 and part != "NA" and index != 1 and float(part) < count_filter:
					line_parts[index] = "NA"
					line_parts[index - 1] = "NA"

			newline = line_parts[0] + "," + line_parts[1]
			for index, part in enumerate(line_parts):
				if(index != 0 and index != 1):
					newline = newline + "," + part

			# Passes the filter if any one of the samples pass both the relative filter and count filter
			if relative_freq_check and count_check:
				table_filtered.write(newline + "\n")
			else:
				table_reads_thrownout.write(newline + "\n")

		line_num = line_num + 1
	table_filtered.close()
	sort_command = "sort -t, -k1,1 -k3,3nr < " + table.split(".csv")[0] + "_filtered.csv > a.tmp && mv a.tmp " + table.split(".csv")[0] + "_filtered.csv"
	subprocess.call(sort_command, shell=True)


# Makes relative frequency plots for each sample by calling syph_visualizer.py.
# Generates plots for both filtered and non-filtered, from the final_data.csvs. 
def bokeh_freq_plot(script_path, cur_dir, new_dir, relative_freq_filter, count_filter):
	for file in os.listdir(cur_dir):
		if file.endswith("fastq"):
			filename = file.split(".fastq")[0]
			print("Constructing relative frequency plots for " + file + "...")
			
			# Grabs sample name
			if file in pacbio_samples_list:
				index = pacbio_samples_list.index(file)
				sample_name = sample_names_list[index] + "_PacBio"
				# PacBio final data will be in this format
				sample_file = filename + ".noprimers.filtered.RAD.nolines.fix_final_data.csv"

			elif file in illumina_samples_list:
				index = illumina_samples_list.index(file)
				sample_name = sample_names_list[index] + "_Illumina"
				# Illumina final data will be in this format
				sample_file = filename + "_final_data.csv"

			# Calls visualizer for non-filtered sample.
			subprocess.call("python " + script_path + "/syph_visualizer.py " + cur_dir + "/" + sample_file + 
			" -t " + sample_name + " -o " + new_dir + " " + svg_flag, shell=True)

			# Filters the current final_data.csv.
			filter_table(relative_freq_filter, count_filter, cur_dir + "/" + sample_file)
			# Recalculates new relative frequencies now that it's filtered.
			subprocess.call("rscript " + script_path + "/recalculate_frequency.R -d " + cwd + " -f " + 
				cur_dir + "/" + sample_file.split(".csv")[0] + "_filtered.csv", shell=True)

			# Calls visualizer for filtered sample.
			subprocess.call("python " + script_path + "/syph_visualizer.py " + cur_dir + "/" + 
				sample_file.split(".csv")[0] + "_filtered.csv -t " + sample_name + "_filtered" + " -o " + new_dir + " " + svg_flag, shell=True)


if __name__ == '__main__': 
	parser = argparse.ArgumentParser(description='tprK project. Currently the program expects '
	 'matched single-end trimmed Illumina and PacBio Q20 reads for the same sample in one big folder.')
	# parser.add_argument('-m', '--metadata_file', required=False,
	# 	help='Specify name of metadata.csv file containing sample name, illumina, and pacbio files. '
	# 	'Otherwise, will look for metadata.csv.')
	parser.add_argument('-f', '--relative_freq_filter', required=False,
		help='Specify by what relative frequency an additional filtered final merged table and visualizations should be sorted at. '
		'By default, this is set to 0.2.')
	parser.add_argument('-c', '--count_filter', required=False,
		help='Specify by what count an additional filtered final merged table and visualizations should be sorted at. '
		'By default, this is set to 5.')
	parser.add_argument('-i', '--illumina_filter', action='store_true', required=False,
		help='Specify if PacBio reads should only include Illumina-supported reads that pass the filters given. '
		'By default, relative freq is set to 0.2 and count is set to 5.')
	parser.add_argument('-pacbio', action='store_true', required = False, help='Write this flag to specify '
		'that there are only PacBio files here. Comparison figures to Illumina will not be created.')
	parser.add_argument('-illumina', action='store_true', required = False, help='Write this flag to specify '
		'that there are only Illumina files here. Comparison figures to PacBio will not be created.')
	parser.add_argument('-svg', action='store_true', required = False, help='Write this flag to specify '
		'.svg files to be created for relative frequency plots for each sample.')

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)

	# Setting variables 
	metadata_file = "metadata.csv"
	if args.relative_freq_filter:
		relative_freq_filter = float(args.relative_freq_filter)
	else:
		relative_freq_filter = 0.2
	if args.count_filter:
		count_filter = int(args.count_filter)
	else:
		count_filter = 5
	if args.pacbio:
		pacbio_flag = "--pacbio"
	else:
		pacbio_flag = ""
	if args.illumina:
		illumina_flag = "--illumina"
	else:
		illumina_flag = ""
	if args.svg:
		svg_flag = "-svg"
	else:
		svg_flag = ""

	script_path = os.path.dirname(os.path.realpath(__file__))
	cwd = os.getcwd()

	sample_names_list = []
	illumina_samples_list = []
	pacbio_samples_list = []

	read_metadata_file = open(metadata_file, "r")
	metadata_file_lines = read_metadata_file.readlines()

	for paired_samples in metadata_file_lines:
		sample, illumina_sample, pacbio_sample = paired_samples.split(",")
		pacbio_sample = pacbio_sample.rstrip()

		# TODO: Run trimmomatic here on illumina_sample
		sample_names_list.append(sample)
		illumina_samples_list.append(illumina_sample)
		pacbio_samples_list.append(pacbio_sample)


	print("Will filter final products for > " + str(relative_freq_filter) + " relative frequency and > " + str(count_filter) + " count.")
	if args.illumina_filter:
		print("Will only include PacBio reads supported by Illumina reads that pass the filter.")

	# Goes from the original Illumina and PacBio reads into the all-important allreads.csv. 
	# Along the way, also makes frequency tables for each sample. 
	print("Running og_files_to_all_reads.py...\n")
	subprocess.call("rscript " + script_path + "/og_files_to_all_reads.R -s " + script_path + " -d " + cwd + " " + 
		pacbio_flag + " " + illumina_flag, shell=True)

	# Creates allreads_filtered.csv and recalculates the relative frequencies.
	filter_table(relative_freq_filter, count_filter, cwd + "/allreads.csv")
	subprocess.call("rscript " + script_path + "/recalculate_frequency.R -d " + cwd + " -f " + cwd + "/allreads_filtered.csv", shell=True)

	# Creates a subfolder in current directory called Figures for the figures to go into.
	figure_output_path = cwd + "/Figures"
	if not os.path.exists("Figures"):
		os.mkdir("Figures")

	# Creates Bokeh frequency plots from the frequency files (all the final_data.csvs).
	# Currently this generates .html plots only, for both filtered and non-filtered for each sample.
	if not args.illumina:
		bokeh_freq_plot(script_path, cwd + "/PacBio_frequencies", figure_output_path, relative_freq_filter, count_filter)
	if not args.pacbio:
		bokeh_freq_plot(script_path, cwd + "/Illumina_frequencies", figure_output_path, relative_freq_filter, count_filter)

	# Generates PacBio vs. Illumina scatterplots for each sample. Compares filtered and non-filtered side by side,
	# as well as automatically generates a zoomed in version from 0-10% relative frequency.
	# Does not occur if running only PacBio or only Illumina files.
	if (not args.pacbio and not args.illumina):
		for num in range(1, len(sample_names_list)):
			print("Generating PacBio vs. Illumina plots for " + sample_names_list[num] + "...")
			subprocess.call("rscript " + script_path + "/PacBio_v_Illumina_plots.R -p " + cwd + 
				" -s " + sample_names_list[num], shell=True)
	else:
		print("-pacbio or -illumina specified. Will not generate PacBio vs. Illumina plots.")

	# Creates a subfolder in the Figures folder for the variable region comparisons to attempt to contain madness.
	comparison_plots_path = cwd + "/Figures/Variable_Region_Comparisons/"
	if not os.path.exists("Figures/Variable_Region_Comparisons"):
		os.mkdir(comparison_plots_path)

	# Creates dot-line plots comparing variable regions between two samples at a time.
	# Lines connect the same tprK region in both samples. Dots indicate that read exists only in that sample.
	# Currently this uses Illumina data.
	if not args.pacbio:
		print("Generating variable region comparison plots...")
		subprocess.call("rscript " + script_path + "/Variable_region_compare.R -p " + cwd, shell=True)

	# Creates a tree based off the filtered cumulative data (alldata_filtered.csv). 
	# Currently no rooting takes place (until somebody figures out how to automate that).
	# Currently this uses PacBio data.
	if not args.illumina:
		print("Generating tree...")
		subprocess.call("rscript " + script_path + "/PacBio2Tree.R -d " + cwd, shell=True)
