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
from math import pi
import subprocess
import os

color_palette=[
	'#5e4fa2', '#3288bd', '#66c2a5', '#abdda4', '#e6f598', '#ffffbf', '#fee08b', '#fdae61', 
	'#f46d43', '#d53e4f', '#9e0142',
	'#003366', '#dec9ab', '#d25757', '#f7d708','#1c29b5', '#b0c997', '#005555', '#f9c3d3', '#0f4c81', 
	'#edc06a', '#bd5915', '#b097c9', '#c2002c', '#808080', 
	'#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3',
	'#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666',
	'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', 
	'#d9d9d9','#bc80bd', '#ccebc5', '#ffed6f',

	'#bc9797', '#a0db8e', '#d69ce1', '#caffcd', '#ffcaf8', '#cafff7', '#f1b4b2', '#030449', 

	'#feff97', '#fd81eb', '#8fb8ff', '#bdffa3', '#ffb5b5', '#0f4c81', '#d7d8ee','#73c7de','#aaaacc',
	'#758eb7', '#f2b6ae', '#f0a890', '#a860a8', '#604890', '#f07878', '#443133', '#ffdecc', '#65bfc1',
	'#ffad19', '#ff9e7d', '#dda97b','#a1b1cc','#d6a562','#c9a48f','#c7906d','#f9733e', 
	'#ffc104', '#624a4c', '#c1b8c9',

			"#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
			"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
			"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
			"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
			"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
			"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
			"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
			"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
			
			"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
			"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
			"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
			"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
			"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
			"#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
			"#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
			"#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
			
			"#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
			"#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
			"#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
			"#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
			"#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
			"#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
			"#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
			"#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",

			"#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
			"#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
			"#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
			"#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
			"#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
			"#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
			"#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
			"#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
			
			"#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
			"#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B",

			'#5e4fa2', '#3288bd', '#66c2a5', '#abdda4', '#e6f598', '#ffffbf', '#fee08b', '#fdae61', 
	'#f46d43', '#d53e4f', '#9e0142',
	'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', 
	'#d9d9d9','#bc80bd', '#ccebc5', '#ffed6f',
	'#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3',
	'#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666',

	'#bc9797', '#a0db8e', '#d69ce1', '#caffcd', '#ffcaf8', '#cafff7', '#f1b4b2', '#030449', 
	'#003366', '#dec9ab', '#d25757', '#f7d708','#1c29b5', '#b0c997', '#005555', '#f9c3d3', '#0f4c81', 
	'#edc06a', '#bd5915', '#b097c9', '#c2002c', '#808080', 
	'#feff97', '#fd81eb', '#8fb8ff', '#bdffa3', '#ffb5b5', '#0f4c81', '#d7d8ee','#73c7de','#aaaacc',
	'#758eb7', '#f2b6ae', '#f0a890', '#a860a8', '#604890', '#f07878', '#443133', '#ffdecc', '#65bfc1',
	'#ffad19', '#ff9e7d', '#dda97b','#a1b1cc','#d6a562','#c9a48f','#c7906d','#f9733e', 
	'#ffc104', '#624a4c', '#c1b8c9',


			"#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
			"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
			"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
			"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
			"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
			"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
			"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
			"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
			
			"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
			"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
			"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
			"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
			"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
			"#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
			"#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
			"#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
			
			"#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
			"#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
			"#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
			"#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
			"#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
			"#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
			"#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
			"#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",

			"#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
			"#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
			"#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
			"#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
			"#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
			"#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
			"#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
			"#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
			
			"#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
			"#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"
    ]

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
			plot_height = 900, plot_width = 1500, min_border_left = 200,
			title = variable_region, toolbar_location = None, sizing_mode = "scale_width")
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

		# Rotate x-axis labels
		fig.xaxis.major_label_orientation = pi/4

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
		mapper = LinearColorMapper(palette=brew_pal, low=2, high=100, low_color = "#ededed")

		hm = figure(
			x_range=sample_list, y_range=sample_reads, title = variable_region, 
			toolbar_location = None, x_axis_location = "above", #background_fill_color = "#d3d3d3",
			#plot_height = (len(sample_reads)*20), 
			plot_height = 350,
			plot_width = (len(sample_list)*25) + (len(max(sample_reads, key = len)) * 8),
			#plot_height = 500,
			#aspect_scale = 1, match_aspect = True,
			#sizing_mode = "scale_both",
			min_border_right = 80,
			min_border_bottom = 80,
			y_axis_location = "left",
			tooltips=[('Read','@y'),('Strain','@x'),('Frequency','@values'+"%")])
		hm.rect(x='x', y='y', width = 1, height = 1, source=source,
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
		hm.background_fill_color = None
		hm.border_fill_color = None
		hm.title.text_font_size = "16pt"

		variable_region_hms.append(hm)

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
		save(grid2)