# -*- coding: iso-8859-1 -*-
# plot_results.py
# DC simulation methods
# Copyright 2009 Giuseppe Venturini

# This file is part of the ahkab simulator.
#
# Ahkab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 of the License.
#
# Ahkab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v2
# along with ahkab.  If not, see <http://www.gnu.org/licenses/>.

"""
This module offers the functions needed to plot the results
of a simulation
"""

import re
import Gnuplot, numpy
import printing, options

def read_data_header(filename):
	fp = open(filename, "r")
	line = fp.readline()
	if line[0] == '#':
		labels = line[1:].upper().split()
	else:
		printing.print_general_error("Unrecognized file: "+filename)
	fp.close()
	return labels

def get_data_label_index(label, filename, labels=None):
	if labels is None:
		labels = read_data_header(filename)
	return labels.index(label.upper())

def read_data(filename, label, labels=None):
	if labels == None:
		labels = read_data_header(filename)
	try:
		index = labels.index(label)
		fp = open(filename, "r")
		data = []
		for line in fp:
			if line.strip()[0] != '#':
				data.append(float(line.split()[index]))
		fp.close()
		data = numpy.array(data)
	except ValueError:
		printing.print_general_error("Unrecognized data set: "+label)
		data =None
	return data

def split_netlist_label(labels_string):
	labels_string = labels_string.strip().upper()
	ret_labels = []
	p = re.compile(r'V\s*\(\s*(\w*)\s*,\s*(\w*)\s*\)', re.IGNORECASE)
	labels_list = p.findall(labels_string)
	for i in range(len(labels_list)):
		l2 = "V"+labels_list[i][0]
		l1 = "V"+labels_list[i][1]
		ret_labels.append((l2,l1))
	p = re.compile(r'V\s*\(\s*(\w*)\s*\)', re.IGNORECASE)
	labels_list = p.findall(labels_string)
	for i in range(len(labels_list)):
		l2 = "V"+labels_list[i]
		l1 = None
		ret_labels.append((l2,l1))
	p = re.compile(r'I\s*\(\s*(\w*)\s*\)', re.IGNORECASE)
	labels_list = p.findall(labels_string)
	for i in range(len(labels_list)):
		l2 = "I("+labels_list[i]+")"
		l1 = None
		ret_labels.append((l2,l1))
	if len(ret_labels) == 0:
		raise Exception, "Unrecognized plot labels: "+ label
	return ret_labels			

def setup_plot(title, xlabel, y2y1_list, analysis):
	g = Gnuplot.Gnuplot()
	g.title(title)
	if xlabel == 'T':
		g.xlabel('t [s]')
	elif xlabel[0] == 'V':
		g.xlabel(xlabel + " [V]")
	elif xlabel[0] == 'I':
		g.xlabel(xlabel + " [A]")
	# here we hope all variables are of the same type
	if y2y1_list[0][0] == 'V':
		g.ylabel("V [V]")
	elif y2y1_list[0][0] == 'I':
		g.ylabel("I [A]")
	return g

def plot_file(title, x, y2y1_list, filename, analysis, outfilename):
	g = setup_plot(title, x, y2y1_list, analysis)
	gfiles = []
	x_index =  str(get_data_label_index(x, filename) + 1)

	plotting_strings = []
	for y2label, y1label in y2y1_list:
		y2_index =  str(get_data_label_index(y2label, filename) + 1)
		if y1label is not None:
			y1_index =  str(get_data_label_index(y1label, filename) + 1)
			plotting_strings.append("($" + y2_index +"-$"+y1_index+")")
		else:
			plotting_strings.append(y2_index)
	
	files_list = []
	for pstring in plotting_strings:
		f = Gnuplot.File(filename, using=x_index+":"+pstring, with=options.plotting_style)
		files_list.append(f)
	
	g.plot(*files_list)

	if outfilename is not None and options.plotting_outtype is not None:
		g.hardcopy(outfilename, terminal=options.plotting_outtype)
	
	if options.plotting_wait_after_plot:
		raw_input('Please press return to continue...\n')
	g.reset()
	del g, files_list

def plot_data(title, x, y2y1_list, filename, analysis, outfilename):
	"""Reads out data from a file and plots it.
	Deprecated, use plot_file, it passes the file directly to gnuplot.
	"""
	
	g = setup_plot(title, x, y2y1_list, analysis)
	gdata = []
	gx = read_data(filename, x)
	for y2label, y1label in y2y1_list:
		if y1label is not None and y1label != '':
			data1 = read_data(filename, y1label)
			ylabel = y2label+"-"+y1label
		else:
			ylabel = y2label
			data1 = 0
		data2 =  read_data(filename, y2label)
		d = Gnuplot.Data(gx, data2-data1, title=ylabel+" ("+analysis+")", with_=options.plotting_style)
		gdata.append(d)
	
	g.plot(*gdata)

	if outfilename is not None and options.plotting_outtype is not None:
		g.hardcopy(outfilename, terminal=options.plotting_outtype)
	
	if options.plotting_wait_after_plot:
		raw_input('Please press return to continue...\n')
	g.reset()
	del data1, data2, g, d, gdata

if __name__ == '__main__':
	filename = 'colpitts_graph.tran'
	labels = read_data_header(filename)
	print labels
	for label in labels[1:]:
                print label
		plot_file("Plot "+label, labels[0], [label], [None], filename, "tran", "plot-"+label+".png")
