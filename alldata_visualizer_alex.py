import numpy as np
import pandas as pd
import argparse
import sys
from bokeh import events
from bokeh.io import save, export_svgs, export_png
from bokeh.resources import CDN
from bokeh.embed import components, file_html
from bokeh.palettes import brewer
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, LabelSet, ColorBar, BasicTicker, PrintfTickFormatter,HoverTool, LinearColorMapper, Slider, CustomJS, Label, WheelZoomTool, ResetTool, Button, TextInput
from bokeh.models.widgets import Panel, Paragraph, Div
from bokeh.layouts import gridplot, column, layout, widgetbox, row
from bokeh.transform import transform
from palette import color_palette
import subprocess
import os

def configure_plot(fig, sample_list, sample_xaxis):
	fig.xaxis.major_label_overrides = dict(zip(sample_xaxis, sample_list))
	fig.xaxis.minor_tick_line_color = None
	fig.title.text_font_size = "28pt"
	fig.xaxis.major_label_text_font_size = '16pt'
	fig.yaxis.major_label_text_font_size = '16pt'
	fig.yaxis.axis_label = "Relative Frequency"
	fig.yaxis.axis_label_text_font_size = '18pt'

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Syph visualizer for all data')
	parser.add_argument('alldata_filtered', help='alldata_filtered.csv data table from tprk_pipeline.py')
	parser.add_argument('metadata', help='metadata.csv, same as used with tprk_pipeline.py')
	parser.add_argument('-svg', action='store_true', help='Use this flag to output graphs in .svg format. '
		'By default, plots will be in .html.')
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	df = pd.read_csv(args.alldata_filtered,index_col=False)
	metadata_file = args.metadata
	source = ColumnDataSource(df)

	# Get list of sample names from metadata file and make xaxis spots for all of them.
	sample_list = []
	sample_xaxis = []
	read_metadata_file = open(metadata_file, "r")
	metadata_file_lines = read_metadata_file.readlines()
	for paired_samples in metadata_file_lines:
		sample, illumina_sample, pacbio_sample = paired_samples.split(",")
		if(sample!="SampleName"):
			sample_list.append(sample)
			sample_xaxis.append(len(sample_list)*10)

	variable_regions_list = ['V1','V2','V3','V4','V5','V6','V7']
	variable_region_figs = []
	variable_region_hms = []
	df_columns = list(df)
	sample_list.sort()	

	# Drop PacBio columns
	for column in df_columns:
		if "PB_" in column or "Count" in column:
			df = df.drop(labels=column, axis=1)
	df_columns = list(df)
	# Replace NaNs with 0s
	df.fillna(0, inplace=True)

	# Initialize palettes
	blues = brewer['Blues'][9][0:7]
	bugn = brewer['BuGn'][9][0:7]
	bupu = brewer['Purples'][8][0:7]
	orrd = brewer['OrRd'][9][0:7]
	gnbu = brewer['GnBu'][9][0:7]
	purd = brewer['PuRd'][9][0:7]
	ylgn = brewer['YlGn'][9][0:7]


	for variable_region in variable_regions_list:
		region_df = df.loc[df.Region==variable_region]
		region_df = region_df.drop(labels="Region", axis=1)
		color_num = 0
		sample_reads = []
		hm_x = []
		hm_y = []
		hm_values = []
		hm_yaxis = []

		fig = figure(x_range = sample_list, y_range = (0,105), tools = "hover", tooltips = "$name",
			plot_height = 1000, plot_width = 1600, title = variable_region, toolbar_location = None)
		data = {'samples': sample_list}
		configure_plot(fig, sample_list, sample_xaxis)
		for index, row in region_df.iterrows():
			read_seq = row[0]
			sample_frequencies = []
			num = 0
			row_parts_all_zero = True
			for i, row_parts in enumerate(row):
				if(i!=0):
					num = float(row_parts)
					if(float(num)!=0):
						row_parts_all_zero = False
						hm_x.append(sample_list[i-1])
						hm_y.append(read_seq)
						hm_values.append(num)
					sample_frequencies.append(num)
			# Only use rows with actual data in it
			if not row_parts_all_zero and read_seq!="":
				sample_reads.append(read_seq)
				data[read_seq] = sample_frequencies				
				color_num = color_num + 1
		fig.vbar_stack(sample_reads, x = 'samples', width = 0.9, source=data,
			fill_color=color_palette[0:len(sample_reads)], fill_alpha = 1, hover_alpha = 0.6,
			hover_color = color_palette[0:len(sample_reads)], line_alpha = 0,
			)
		variable_region_figs.append(fig)

		df2 = pd.DataFrame()
		df2['x']=hm_x
		df2['y']=hm_y
		df2['values']=hm_values
		source = ColumnDataSource(df2)
		if(variable_region=="V1"):
			brew_pal = blues
		elif(variable_region=="V2"):
			brew_pal = bugn 
		elif(variable_region=="V3"):
			brew_pal = bupu
		elif(variable_region=="V4"):
			brew_pal = orrd 
		elif(variable_region=="V5"):
			brew_pal = gnbu 
		elif(variable_region=="V6"):
			brew_pal = ylgn 
		else:
			brew_pal = purd 
		# Reverse the colors so darker colors at max
		brew_pal = brew_pal[::-1]
		mapper = LinearColorMapper(palette=brew_pal, low=2, high=100, low_color = "#ffffff")

		hm = figure(x_range=sample_list, y_range=sample_reads, title = variable_region, 
			toolbar_location = None, x_axis_location = "above", #background_fill_color = "#d3d3d3",
			plot_height = (len(sample_reads)*20), 
			plot_width = (len(sample_list)*25) + (len(max(sample_reads, key = len)) * 8),
			min_border_right = 80,
			y_axis_location = "left",
			tooltips=[('Read','@y'),('Strain','@x'),('Frequency','@values')])
		hm.rect(x='x', y='y', width = 1, height =1, source=source,
			line_color=None, fill_color=transform('values', mapper), hover_line_color = 'black',
			hover_color = transform('values',mapper))
		color_bar = ColorBar(color_mapper = mapper, major_label_text_font_size="5pt",
			ticker=BasicTicker(desired_num_ticks=len(purd)),
			formatter=PrintfTickFormatter(format="%d%%"),
			label_standoff=6, border_line_color=None, location=(0, 0))
		hm.add_layout(color_bar, 'right')
		hm.xaxis.major_label_orientation = 1.0
		hm.grid.grid_line_color = None
		hm.axis.axis_line_color = None
		hm.axis.major_tick_line_color = None
		variable_region_hms.append(hm)

		hm.background_fill_color = None
		hm.border_fill_color = None


	grid = gridplot(variable_region_figs, ncols=1)
	output_file("all_variable_regions.html")
	save(grid)


	if args.svg:
		for myfig in variable_region_hms:
			print(myfig.title.text, " ", myfig.plot_height)
			myfig.output_backend = "svg"
			output_filename = myfig.title.text + ".svg"
			output_file(output_filename)
			export_svgs(myfig, filename=output_filename)

	else:
		grid2 = gridplot(variable_region_hms,ncols=1)
		output_file("all_heatmap.html", title = "Heatmap")
		save(variable_region_hms, sizing_mode = "scale_width")