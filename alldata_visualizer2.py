import numpy as np
import pandas as pd
import argparse
import sys
from bokeh import events
from bokeh.io import save, export_svgs
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
	fig.title.text_font_size = "20pt"
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
	blues = brewer['Blues'][9]
	bugn = brewer['BuGn'][9]
	bupu = brewer['Purples'][9]
	orrd = brewer['OrRd'][9]
	gnbu = brewer['GnBu'][9]
	purd = brewer['PuRd'][9]
	ylgn = brewer['YlGn'][9]
	rdpu = brewer['RdPu'][9]
	ylbu = brewer['YlGnBu'][9]
	ylbr = brewer['YlOrBr'][9]

	brew_pal_list = blues, bugn, bupu, orrd, gnbu, ylgn, purd,rdpu,ylbu,ylbr
	xaxis_locations = []

	for variable_region in variable_regions_list:
		xaxis_locations = []
		region_df = df.loc[df.Region==variable_region]
		region_df = region_df.drop(labels="Region", axis=1)

		fig = figure(x_range = (1,(len(sample_list)*10)+9), y_range = (0,105),
			plot_height = 1000, plot_width = 1600, title = variable_region, toolbar_location = None)

		for index, sample in enumerate(sample_list):
			sample_col = list(region_df.iloc[:,index+1])
			# Removes 0s
			sample_col = [x for x in sample_col if float(x)!= 0]
			sample_col = sorted(sample_col, reverse=True)
			palette = brew_pal_list[index%10]

			sample_bottom = 0
			read_count = 0
			for read in sample_col:
				bar = fig.vbar(x=10*(index+1), width = 6, bottom = sample_bottom, top = [sample_bottom + float(read)], color=palette[read_count%9])
				read_count = read_count + 1
				sample_bottom = sample_bottom + float(read)
			xaxis_locations.append(10*(index+1))

		fig.xaxis.ticker = xaxis_locations
		fig.xaxis.major_label_overrides = dict(zip(xaxis_locations, sample_list))
		fig.xaxis.minor_tick_line_color = None
		configure_plot(fig, sample_list,xaxis_locations)

		variable_region_figs.append(fig)

	grid = gridplot(variable_region_figs, ncols=1)
	output_file("all_variable_regions2.html")
	save(grid)