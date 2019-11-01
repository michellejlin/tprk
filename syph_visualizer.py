import numpy as np
import pandas as pd
import argparse
import sys
from bokeh.palettes import brewer
from bokeh import events
from bokeh.io import export_png, save, export_svgs
from bokeh.resources import CDN
from bokeh.embed import components, file_html, autoload_static
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, Jitter, HoverTool, Slider, CustomJS, Label, WheelZoomTool, ResetTool, Button, TextInput
from bokeh.transform import jitter, factor_cmap
from bokeh.models.widgets import Panel, Tabs, Paragraph, Div, CheckboxGroup
from bokeh.layouts import column, layout, widgetbox, row
import subprocess
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Syph visualizer for single isolates')
	parser.add_argument('final_data', help='data table from syph.py, in csv format.')
	parser.add_argument('-t', '--title', help='name of sample to use as plot title.')
	parser.add_argument('-o', '--output', help='output folder to export plots.')
	parser.add_argument('-svg', action='store_true', help='Use this flag to output graphs in .svg format. '
		'By default, plots will be in .html.')
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	# Setting variables
	# Reads data from syph.py output.
	strain = args.title
	table = pd.read_csv(args.final_data,index_col=False)
	source = ColumnDataSource(table)

	#TOOLTIPS = [
	#	("Read", "@Read"),
	#	("Count", "@Count"),
	#	("Relative Frequency", "@RelativeFreq"+"%"),
	#]

	plot_title = strain
	output_path = args.output

	if args.svg:
		fig1 = figure(x_range = (0,80), y_range = (0,105), plot_height = 330, plot_width = 540, title = strain, toolbar_location = None)
	else:
		fig1 = figure(x_range = (0,80), y_range = (0,105), plot_height = 1000, plot_width = 1600, title = strain, toolbar_location = None)
	
	#TODO: Fix 0 on x-axis
	fig1.xaxis.major_label_overrides = dict(zip([10,20,30,40,50,60,70,80], ["V1", "V2", "V3", "V4", "V5", "V6", "V7",""]))
	fig1.xaxis.minor_tick_line_color = None
	fig1.title.text_font_size = "18pt"
	fig1.xaxis.major_label_text_font_size = '16pt'
	fig1.yaxis.major_label_text_font_size = '16pt'
	fig1.yaxis.axis_label = "Relative Frequency"
	fig1.yaxis.axis_label_text_font_size = '18pt'

	blues = brewer['Blues'][9]
	bugn = brewer['BuGn'][9]
	bupu = brewer['Purples'][8]
	orrd = brewer['OrRd'][9]
	gnbu = brewer['GnBu'][9]
	purd = brewer['PuRd'][9]
	ylgn = brewer['YlGn'][9]
	v1_bottom = 0
	v2_bottom = 0
	v3_bottom = 0
	v4_bottom = 0
	v5_bottom = 0
	v6_bottom = 0
	v7_bottom = 0
	region_totals = [0,0,0,0,0,0,0]
	region_counts = [0,0,0,0,0,0,0]

	for line in open(args.final_data):
		region, read, relativefreq, count = line.split(',')
		if (region == "V1"):
			bar = fig1.vbar(x=10, width = 6, bottom = v1_bottom, top = [v1_bottom + float(relativefreq)], 
				color=blues[region_counts[0]%9])
			hover = HoverTool(renderers = [bar], toggleable = False,
				tooltips=[("Read", read), ("Count", count), ("Relative Frequency", relativefreq)])
			fig1.add_tools(hover)
			v1_bottom = v1_bottom + float(relativefreq)
			region_totals[0] = region_totals[0] + int(count)
			region_counts[0] = region_counts[0] + 1

		elif (region == "V2"):
			bar = fig1.vbar(x=20, width = 6, bottom = v2_bottom, top = [v2_bottom + float(relativefreq)], 
				color=bugn[region_counts[1]%9])
			hover = HoverTool(renderers = [bar], toggleable = False,
				tooltips=[("Read", read), ("Count", count), ("Relative Frequency", relativefreq)])
			fig1.add_tools(hover)
			v2_bottom = v2_bottom + float(relativefreq)
			region_totals[1] = region_totals[1] + int(count)
			region_counts[1] = region_counts[1] + 1
		elif (region == "V3"):
			bar = fig1.vbar(x=30, width = 6, bottom = v3_bottom, top = [v3_bottom + float(relativefreq)], 
				color=bupu[region_counts[2]%8])
			hover = HoverTool(renderers = [bar], toggleable = False,
				tooltips=[("Read", read), ("Count", count), ("Relative Frequency", relativefreq)])
			fig1.add_tools(hover)
			v3_bottom = v3_bottom + float(relativefreq)
			region_totals[2] = region_totals[2] + int(count)
			region_counts[2] = region_counts[2] + 1
		elif (region == "V4"):
			bar = fig1.vbar(x=40, width = 6, bottom = v4_bottom, top = [v4_bottom + float(relativefreq)], 
				color=orrd[region_counts[3]%9])
			hover = HoverTool(renderers = [bar], 
				tooltips=[("Read", read), ("Count", count), ("Relative Frequency", relativefreq)])
			fig1.add_tools(hover)
			v4_bottom = v4_bottom + float(relativefreq)
			region_totals[3] = region_totals[3] + int(count)
			region_counts[3] = region_counts[3] + 1
		elif (region == "V5"):
			bar = fig1.vbar(x=50, width = 6, bottom = v5_bottom, top = [v5_bottom + float(relativefreq)], 
				color=gnbu[region_counts[4]%9])
			hover = HoverTool(renderers = [bar], toggleable = False,
				tooltips=[("Read", read), ("Count", count), ("Relative Frequency", relativefreq)])
			fig1.add_tools(hover)
			v5_bottom = v5_bottom + float(relativefreq)
			region_totals[4] = region_totals[4] + int(count)
			region_counts[4] = region_counts[4] + 1
		elif (region == "V6"):
			bar = fig1.vbar(x=60, width = 6, bottom = v6_bottom, top = [v6_bottom + float(relativefreq)], 
				color=ylgn[region_counts[5]%9])
			hover = HoverTool(renderers = [bar], toggleable = False,
				tooltips=[("Read", read), ("Count", count), ("Relative Frequency", relativefreq)])
			fig1.add_tools(hover)
			v6_bottom = v6_bottom + float(relativefreq)
			region_totals[5] = region_totals[5] + int(count)
			region_counts[5] = region_counts[5] + 1
		elif (region == "V7"):
			bar = fig1.vbar(x=70, width = 6, bottom = v7_bottom, top = [v7_bottom + float(relativefreq)], 
				color=purd[region_counts[6]%9])
			hover = HoverTool(renderers = [bar], toggleable = False,
				tooltips=[("Read", read), ("Count", count), ("Relative Frequency", relativefreq)])
			fig1.add_tools(hover)
			v7_bottom = v7_bottom + float(relativefreq)
			region_totals[6] = region_totals[6] + int(count)
			region_counts[6] = region_counts[6] + 1

	for index, total in enumerate(region_totals):
		label = Label(x = (index + 1) * 10, y = 101, text = str(region_counts[index]) + ", " + str(total), border_line_color = None, text_align = 'center', text_font_size = '11pt')
		fig1.add_layout(label)

 
	# fig2 = figure(plot_width = 1600, plot_height = 1000, title = "Fig 2", y_range = (0,105), x_range = table.Region.unique(), tooltips=TOOLTIPS)
	# circle = fig2.circle(x=jitter('Region', width = 0.3, range=fig2.x_range), y='RelativeFreq', size = 15, alpha = 0.8, source=source, color='cornflowerblue')
	# fig2.title.text_font_size = "18pt"
	# fig2.xaxis.major_label_text_font_size = '12pt'
	# fig2.yaxis.major_label_text_font_size = '12pt'
	# fig2.yaxis.axis_label = "Relative Frequency"

	if args.svg:
		fig1.output_backend = "svg"
		output_filename = (output_path + "/" + strain + "_RelativeFreqPlot.svg")
		output_file(output_filename)
		export_svgs(fig1, filename=output_filename)
	else:
		output_filename = (output_path + "/" + strain + "_RelativeFreqPlot.html")
		output_file(output_filename, title=strain)

		save(fig1)

	# subprocess.call("mv " + output_filename + " " + args.output + "/" + output_filename, shell=True)
