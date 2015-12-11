#!/usr/local/bin/python3

# Michael Matschiner, 2014-09-11
# michaelmatschiner@mac.com

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
	print('Python 3 is needed to run this script!')
	sys.exit(0)
import argparse
import textwrap
import re
import math

# Parse the command line arguments.
parser = argparse.ArgumentParser(
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description=textwrap.dedent('''\
	  %(prog)s
	-----------------------------------------
	  Extracts the haplotype genealogy or
	  individual statistics from Fitchi
	  HTML files.
	'''))
parser.add_argument(
	'-v', '--version',
	action='version',
	version='%(prog)s 0.97'
	)
parser.add_argument(
    '-e',
    nargs=1,
    type=str,
    default=[None],
    dest='extract',
    help="Section or statistic to be extracted from Fitchi HTML file. \
    	Possible arguments are: svg, svg_simple, svg_bw, \
    	svg_simple_bw (the haplotype genealogy SVG code, see manual for \
    	details), legend (the figure legend), prop_var (the proportion of \
    	variable sites in the alignments), tot_var (the total number of \
    	variable sites), fst, dxy, and df (each for the first pairwise \
    	comparison), as well as the gsi of the first (gsi1) and second \
		(gsi2) population."
    	)
parser.add_argument(
	'infile',
	nargs='?',
	type=argparse.FileType('r'),
	default='-',
	help='The input file name.'
	)
parser.add_argument(
	'outfile',
	nargs='?',
	type=argparse.FileType('w'),
	default=sys.stdout,
	help='The output file name.'
	)
args = parser.parse_args()
extract = args.extract[0]
infile = args.infile
outfile = args.outfile
if infile.isatty():
	print('No input file specified, and no input piped through stdin!')
	sys.exit(0)
instring = infile.read()

if extract in ['svg', 'svg_simple', 'svg_bw', 'svg_simple_bw']:

	# Prepare an array for the SVG lines.
	svg_lines = ['<?xml version="1.0" encoding="utf-8"?>']
	svg_lines.append('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">')

	# parse the population lines and the svg lines.
	lines = instring.split('\n')
	info_start_pattern = re.compile("<!-- The haplotype genealogy graph: info start -->")
	info_end_pattern = re.compile("<!-- The haplotype genealogy graph: info end -->")
	svg_start_pattern = re.compile("<!-- The haplotype genealogy graph: SVG string start -->")
	svg_end_pattern = re.compile("<!-- The haplotype genealogy graph: SVG string end -->")
	in_svg_part = False
	in_info_part = False
	pop_lines = []
	for line in lines:
		svg_hit_start = re.search(svg_start_pattern, line)
		svg_hit_end = re.search(svg_end_pattern, line)
		if svg_hit_start != None:
			in_svg_part = True
		elif svg_hit_end != None:
			in_svg_part = False
		elif in_svg_part == True:
			svg_lines.append(line[12:])
		info_hit_start = re.search(info_start_pattern, line)
		info_hit_end = re.search(info_end_pattern, line)
		if info_hit_start != None:
			in_info_part = True
		elif info_hit_end != None:
			in_info_part = False
		elif in_info_part == True:
			pop_lines.append(line[12:])
	pops = []
	colors = []
	for pop_line in pop_lines[1:-1]:
		pop_line_ary = pop_line.split(': #')
		pops.append(pop_line_ary[0])
		colors.append(pop_line_ary[1])

	# Get the size of the SVG graph without legend.
	dim_x_pattern = re.compile('width=\".+?\"')
	dim_y_pattern = re.compile('height=\".+?\"')
	dim_x_hit = re.search(dim_x_pattern, svg_lines[2])
	dim_y_hit = re.search(dim_y_pattern, svg_lines[2])
	dim_x = int(dim_x_hit.group(0).split('=')[1][1:-1])
	dim_y = int(dim_y_hit.group(0).split('=')[1][1:-1])

	# Define parameters for the legend.
	font_unit = 10
	margin = 10
	stroke_width = 1
	font = 'helvetica'

	# From the population lines, calculate the size of the legend.
	max_pop_length = 0
	for pop in pops:
		if len(pop) > max_pop_length:
			max_pop_length = len(pop)
	legend_size_x = math.ceil(max_pop_length * font_unit + 2 * margin)

	# Change the third line of the svg text (which has the dimensions)
	# to accomodate the legend on the right hand side of the graph.
	extended_dim_x = dim_x + legend_size_x
	svg_third_line = '<svg xmlns="http://www.w3.org/2000/svg" width="' + str(extended_dim_x)
	svg_third_line += '" height="' + str(dim_y)
	svg_third_line += '" viewBox="0 0 ' + str(extended_dim_x) + ' ' + str(dim_y) + '" xmlns:xlink="htttp://www.w3.org/1999/xlink">'
	svg_lines[2] = svg_third_line

	# Prepare the SVG legend list.
	svg_legend_string = ''
	svg_legend_string += '  <!--Legend-->\n'
	svg_legend_string += '  <line x1="' + str(dim_x) + '" y1="0" x2="' + str(dim_x) + '" y2="' + str(dim_y) + '"/>\n'
	for x in range(0,len(pops)):
		svg_legend_string += '  <circle fill="#' + colors[x] + '" cx="'
		svg_legend_string += str(dim_x + margin + 0.5 * (font_unit - stroke_width))
		svg_legend_string += '" cy="'
		svg_legend_string += str(margin + 1.3 * font_unit * x)
		svg_legend_string += '" r="'
		svg_legend_string += str(0.5 * (font_unit- stroke_width))
		svg_legend_string += '"/>\n'
		svg_legend_string += '<text x="'
		svg_legend_string += str(dim_x + margin + 0.5 * (font_unit - stroke_width) + font_unit)
		svg_legend_string += '" y="'
		svg_legend_string += str(margin + 0.35 * font_unit + 1.3 * font_unit * x)
		svg_legend_string += '" font-family="'
		svg_legend_string += font
		svg_legend_string += '" font-size="'
		svg_legend_string += str(font_unit)
		svg_legend_string += 'px">'
		svg_legend_string += pops[x]
		svg_legend_string += '</text>\n'
	svg_legend_string += '\n'
	if extract in ['svg_bw', 'svg_simple_bw']:
		svg_legend_string = svg_legend_string.replace('859900','000000')
		svg_legend_string = svg_legend_string.replace('b58900','dfdfdf')
		svg_legend_string = svg_legend_string.replace('2aa198','6f6f6f')
		svg_legend_string = svg_legend_string.replace('268bd2','2f2f2f')
		svg_legend_string = svg_legend_string.replace('6c71c4','afafaf')
		svg_legend_string = svg_legend_string.replace('d33682','4f4f4f')
		svg_legend_string = svg_legend_string.replace('dc322f','8f8f8f')
		svg_legend_string = svg_legend_string.replace('cb4b16','0f0f0f')
		svg_legend_string = svg_legend_string.replace('002b36','cfcfcf')
		svg_legend_string = svg_legend_string.replace('586e75','1f1f1f')
		svg_legend_string = svg_legend_string.replace('839496','3f3f3f')
		svg_legend_string = svg_legend_string.replace('a2aca8','5f5f5f')
		svg_legend_string = svg_legend_string.replace('c6c6bc','7f7f7f')
		svg_legend_string = svg_legend_string.replace('657b83','9f9f9f')
		svg_legend_string = svg_legend_string.replace('b5bab3','bfbfbf')

	# Turn all lines into a string, making modifications to color
	# and gradients if needed.
	include_svg_line = True
	svg_string = ''
	for svg_line in svg_lines[1:-1]:
		if extract in ['svg_bw', 'svg_simple_bw']:
			svg_legend_string = svg_legend_string.replace('859900','000000')
			svg_legend_string = svg_legend_string.replace('b58900','dfdfdf')
			svg_legend_string = svg_legend_string.replace('2aa198','6f6f6f')
			svg_legend_string = svg_legend_string.replace('268bd2','2f2f2f')
			svg_legend_string = svg_legend_string.replace('6c71c4','afafaf')
			svg_legend_string = svg_legend_string.replace('d33682','4f4f4f')
			svg_legend_string = svg_legend_string.replace('dc322f','8f8f8f')
			svg_legend_string = svg_legend_string.replace('cb4b16','0f0f0f')
			svg_legend_string = svg_legend_string.replace('002b36','cfcfcf')
			svg_legend_string = svg_legend_string.replace('586e75','1f1f1f')
			svg_legend_string = svg_legend_string.replace('839496','3f3f3f')
			svg_legend_string = svg_legend_string.replace('a2aca8','5f5f5f')
			svg_legend_string = svg_legend_string.replace('c6c6bc','7f7f7f')
			svg_legend_string = svg_legend_string.replace('657b83','9f9f9f')
			svg_legend_string = svg_legend_string.replace('b5bab3','bfbfbf')
		if extract in ['svg_simple', 'svg_simple_bw']:
			if '<radialGradient id="radgrad"' in svg_line:
				include_svg_line = False
			elif '</defs>' in svg_line:
				include_svg_line = True
				svg_string += svg_line + '\n'
			elif '<!--Gradients-->' in svg_line:
				include_svg_line = False
			elif '<!--Node labels-->' in svg_line:
				include_svg_line = True
				svg_string += svg_line + '\n'
			elif include_svg_line == True:
				svg_string += svg_line + '\n'
		else:
			svg_string += svg_line + '\n'
	svg_string += svg_legend_string
	svg_string += '</svg>\n'
	outfile.write(svg_string)

# Find the figure legend.
elif extract == 'legend':
	lines = instring.split('\n')
	legend_start_pattern = re.compile("<!-- Legend: string start -->")
	legend_end_pattern = re.compile("<!-- Legend: string end -->")
	in_legend_part = False
	legend_string = ""
	for line in lines:
		legend_hit_start = re.search(legend_start_pattern, line)
		legend_hit_end = re.search(legend_end_pattern, line)
		if legend_hit_start != None:
			in_legend_part = True
		elif legend_hit_end != None:
			in_legend_part = False
		elif in_legend_part == True:
			legend_string += line.strip()
	legend_string = legend_string.replace(" class=\"spaceOver\"","")
	legend_string = legend_string.replace(" class=\"spaceUnder\"","")
	legend_string = legend_string.replace(" width=\"160\"","")
	legend_string = legend_string.replace(" style=\"font-weight:bold; font-family:Courier\"","")
	legend_string = legend_string.replace(" style=\"font-family:Courier\"","")
	legend_string = legend_string.replace(" valign=\"top\"","")
	legend_string = legend_string.replace("<br>","")
	legend_string = legend_string.replace("<div style=\"width: 640px; overflow: auto;\">","")
	legend_string = legend_string.replace("</div>","")
	legend_ary = legend_string.split('</tr><tr>')
	legend_table_string = ""
	for legend_table_row in legend_ary:
		legend_table_row = legend_table_row.replace("<tr>","").replace("</tr>","")
		legend_table_string += legend_table_row.replace("</td><td>","\t").replace("<td>","").replace("</td>","") + "\n"
	outfile.write(legend_table_string)

# Find the total number of variable sites.
elif extract == 'tot_var':
	lines = instring.split('\n')
	tot_var_start_pattern = re.compile("<!-- Total variable: string start -->")
	tot_var_end_pattern = re.compile("<!-- Total variable: string end -->")
	in_tot_var_part = False
	tot_var = None
	for line in lines:
		tot_var_hit_start = re.search(tot_var_start_pattern, line)
		tot_var_hit_end = re.search(tot_var_end_pattern, line)
		if tot_var_hit_start != None:
			in_tot_var_part = True
		elif tot_var_hit_end != None:
			in_tot_var_part = False
		elif in_tot_var_part == True:
			tot_var_pattern = re.compile(">.+<")
			tot_var_string = re.search(tot_var_pattern, line)
			tot_var = int(tot_var_string.group(0)[1:-1])
	outfile.write(str(tot_var) + '\n')

# Find the proportion of variable sites.
elif extract == 'prop_var':
	lines = instring.split('\n')
	prop_var_start_pattern = re.compile("<!-- Proportion variable: string start -->")
	prop_var_end_pattern = re.compile("<!-- Proportion variable: string end -->")
	in_prop_var_part = False
	prop_var = None
	for line in lines:
		prop_var_hit_start = re.search(prop_var_start_pattern, line)
		prop_var_hit_end = re.search(prop_var_end_pattern, line)
		if prop_var_hit_start != None:
			in_prop_var_part = True
		elif prop_var_hit_end != None:
			in_prop_var_part = False
		elif in_prop_var_part == True:
			prop_var_pattern = re.compile(">.+<")
			prop_var_string = re.search(prop_var_pattern, line)
			prop_var = float(prop_var_string.group(0)[1:-1])
	outfile.write(str(prop_var) + '\n')

# Find the probability pi that two randomly drawn individuals have different alleles.
elif extract == 'pi':
	lines = instring.split('\n')
	pi_start_pattern = re.compile("<!-- Pi: string start -->")
	pi_end_pattern = re.compile("<!-- Pi: string end -->")
	in_pi_part = False
	pi = None
	for line in lines:
		pi_hit_start = re.search(pi_start_pattern, line)
		pi_hit_end = re.search(pi_end_pattern, line)
		if pi_hit_start != None:
			in_pi_part = True
		elif pi_hit_end != None:
			in_pi_part = False
		elif in_pi_part == True:
			pi_pattern = re.compile(">.+<")
			pi_string = re.search(pi_pattern, line)
			pi = float(pi_string.group(0)[1:-1])
	outfile.write(str(pi) + '\n')

# Find the first pairwise Fst value.
elif extract == 'fst':
	lines = instring.split('\n')
	fst_start_pattern = re.compile("<!-- First f_st: string start -->")
	fst_end_pattern = re.compile("<!-- First f_st: string end -->")
	in_fst_part = False
	fst = None
	for line in lines:
		fst_hit_start = re.search(fst_start_pattern, line)
		fst_hit_end = re.search(fst_end_pattern, line)
		if fst_hit_start != None:
			in_fst_part = True
		elif fst_hit_end != None:
			in_fst_part = False
		elif in_fst_part == True:
			fst_pattern = re.compile(">.+<")
			fst_string = re.search(fst_pattern, line)
			if fst_string.group(0)[1:-1] == 'NA':
				fst = 'NA'
			else:
				fst = float(fst_string.group(0)[1:-1])
	outfile.write(str(fst) + '\n')

# Find the first d_xy value (see Ruegg et al. 2014).
elif extract == 'dxy':
	lines = instring.split('\n')
	dxy_start_pattern = re.compile("<!-- First d_xy: string start -->")
	dxy_end_pattern = re.compile("<!-- First d_xy: string end -->")
	in_dxy_part = False
	dxy = None
	for line in lines:
		dxy_hit_start = re.search(dxy_start_pattern, line)
		dxy_hit_end = re.search(dxy_end_pattern, line)
		if dxy_hit_start != None:
			in_dxy_part = True
		elif dxy_hit_end != None:
			in_dxy_part = False
		elif in_dxy_part == True:
			dxy_pattern = re.compile(">.+<")
			dxy_string = re.search(dxy_pattern, line)
			if dxy_string.group(0)[1:-1] == 'NA':
				dxy = 'NA'
			else:
				dxy = float(dxy_string.group(0)[1:-1])
	outfile.write(str(dxy) + '\n')

# Find the first d_f value (see Ruegg et al. 2014).
elif extract == 'df':
	lines = instring.split('\n')
	df_start_pattern = re.compile("<!-- First d_f: string start -->")
	df_end_pattern = re.compile("<!-- First d_f: string end -->")
	in_df_part = False
	df = None
	for line in lines:
		df_hit_start = re.search(df_start_pattern, line)
		df_hit_end = re.search(df_end_pattern, line)
		if df_hit_start != None:
			in_df_part = True
		elif df_hit_end != None:
			in_df_part = False
		elif in_df_part == True:
			df_pattern = re.compile(">.+<")
			df_string = re.search(df_pattern, line)
			if df_string.group(0)[1:-1] == 'NA':
				df = 'NA'
			else:
				df = float(df_string.group(0)[1:-1])
	outfile.write(str(df) + '\n')

# Find the first gsi value.
elif extract == 'gsi1':
	lines = instring.split('\n')
	gsi1_start_pattern = re.compile("<!-- First gsi: string start -->")
	gsi1_end_pattern = re.compile("<!-- First gsi: string end -->")
	in_gsi1_part = False
	gsi1 = None
	for line in lines:
		gsi1_hit_start = re.search(gsi1_start_pattern, line)
		gsi1_hit_end = re.search(gsi1_end_pattern, line)
		if gsi1_hit_start != None:
			in_gsi1_part = True
		elif gsi1_hit_end != None:
			in_gsi1_part = False
		elif in_gsi1_part == True:
			gsi1_pattern = re.compile(">.+<")
			gsi1_string = re.search(gsi1_pattern, line)
			if gsi1_string.group(0)[1:-1] == 'NA':
				gsi1 = 'NA'
			else:
				gsi1 = float(gsi1_string.group(0)[1:-1])
	outfile.write(str(gsi1) + '\n')

elif extract == 'gsi2':
	lines = instring.split('\n')
	gsi2_start_pattern = re.compile("<!-- Second gsi: string start -->")
	gsi2_end_pattern = re.compile("<!-- Second gsi: string end -->")
	in_gsi2_part = False
	gsi2 = None
	for line in lines:
		gsi2_hit_start = re.search(gsi2_start_pattern, line)
		gsi2_hit_end = re.search(gsi2_end_pattern, line)
		if gsi2_hit_start != None:
			in_gsi2_part = True
		elif gsi2_hit_end != None:
			in_gsi2_part = False
		elif in_gsi2_part == True:
			gsi2_pattern = re.compile(">.+<")
			gsi2_string = re.search(gsi2_pattern, line)
			if gsi2_string.group(0)[1:-1] == 'NA':
				gsi2 = 'NA'
			else:
				gsi2 = float(gsi2_string.group(0)[1:-1])
	outfile.write(str(gsi2) + '\n')

elif extract == None:
	print("Nothing to be extracted from Fitchi HTML file (use '-h' to see available options).")

