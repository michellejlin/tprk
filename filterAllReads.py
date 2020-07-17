import subprocess
import argparse
import sys
import os 
import regex

# Filters allreads.csv and outputs allreads_filtered.csv based on flags -f and -c.
# (Default relative freq is 0.2 and count is 5.) Thrown away reads are put into allreads_thrownout.csv.
def filter_table(relative_freq_filter, count_filter, table, is_heatmap):
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
			if(is_heatmap):
				print("Generating allreads_filtered_heatmap.csv file...")
				table_filtered = open(table.split(".csv")[0] + "_filtered_heatmap.csv", "w+")
			else:
				table_filtered = open(table.split(".csv")[0] + "_filtered.csv", "w+")
			table_filtered.write(line + "\n")
		# Filters lines based on specified values. Checks that each sample passes both the relative frequency filter
		# and the count filter. Writes passed lines to alldata_filtered.csv, and thrown away lines to alldata_thrownout.csv.
		else:
			line_parts = line.split(",")
			relative_freq_check = False
			count_check = False

			newline = line_parts[0] + "," + line_parts[1]
			newline_hm = line_parts[0] + "," + line_parts[1]

			# print(line_parts[1])

			for index, part in enumerate(line_parts):
				if (index != 0 and index != 1):
					# For non-heatmap, we don't want to include parts that don't pass filters, even if there are other samples with the same read that do.
					if (is_heatmap == False):
						if (index % 2 == 1):
							# print(part)
							# print(line_parts[index-1])
							if (part != "NA" and (float(part) >= count_filter) and (float(line_parts[index - 1]) >= relative_freq_filter)):
								newline = newline + "," + line_parts[index - 1] + "," + part
								relative_freq_check = True
								count_check = True
								# print(index)
								# print(newline)
							else:
								newline = newline + ",NA,NA"
								# print(newline)
					# Heatmap only! Check if there are samples with the same read name that do have the read, in which case we keep the whole line.
					# Even indexes should be relative frequencies
					else:					
						if (part != "NA" and (index % 2) == 0) and float(part) >= relative_freq_filter:
							relative_freq_check = True
						# Odd indexes should be counts
						if (part != "NA" and (index % 2 == 1)) and float(part) >= count_filter:
							count_check = True
						# For allreads_filtered.csv, we want to change anything under the count filter to NA.
						if (part != "NA" and ("allreads" in table) and (index % 2 == 1) and float(part) < count_filter):
							line_parts[index] = "NA"
							line_parts[index - 1] = "NA"

			# Create the whole line again for the heatmap.
			for index, part in enumerate(line_parts):
				if(index != 0 and index != 1):
					newline_hm = newline_hm + "," + part

			# Passes the filter if any one of the samples pass both the relative filter and count filter, only for heatmap
			if relative_freq_check and count_check:
				if is_heatmap:
					table_filtered.write(newline_hm + "\n")
				else:
					table_filtered.write(newline + "\n")
			elif not is_heatmap:
				table_reads_thrownout.write(newline + "\n")

		line_num = line_num + 1
	table_filtered.close()

	if is_heatmap:
		sort_command = "sort -t, -k1,1 -k3,3nr < " + table.split(".csv")[0] + "_filtered_heatmap.csv > a.tmp && mv a.tmp " + table.split(".csv")[0] + "_filtered_heatmap.csv"
	else:
		sort_command = "sort -t, -k1,1 -k3,3nr < " + table.split(".csv")[0] + "_filtered.csv > a.tmp && mv a.tmp " + table.split(".csv")[0] + "_filtered.csv"
	subprocess.call(sort_command, shell=True)

if __name__ == '__main__': 
	parser = argparse.ArgumentParser(description='tprK project. Currently the program expects '
	 'matched single-end trimmed Illumina and PacBio Q20 reads for the same sample in one big folder.')
	parser.add_argument('-a', '--allreads_file', required=False,
		help='Specify name of allreads.csv')
	parser.add_argument('-f', '--relative_freq_filter', required=False,
		help='Specify by what relative frequency an additional filtered final merged table and visualizations should be sorted at. '
		'By default, this is set to 0.2.')
	parser.add_argument('-c', '--count_filter', required=False,
		help='Specify by what count an additional filtered final merged table and visualizations should be sorted at. '
		'By default, this is set to 5.')
	parser.add_argument('-is_heatmap', action='store_true', required = False, help='')

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)
    
	relative_freq_filter = float(args.relative_freq_filter)
	count_filter = float(args.count_filter)

	filter_table(relative_freq_filter, count_filter, args.allreads_file, args.is_heatmap)