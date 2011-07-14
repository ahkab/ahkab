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
import numpy, pylab
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
	label = label.upper()
 	if labels == None:
 		labels = read_data_header(filename)
	else:
		labels = map(str.upper, labels)
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
		data = None
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
		raise Exception, "Unrecognized plot labels: "+ labels_string
	return ret_labels			

def setup_plot(fig, title, xlabel, y2y1_list, analysis, log=False):
	pylab.title(title)
	ax = pylab.gca()
	if xlabel == 'T':
		pylab.xlabel('t [s]')
	elif xlabel[0] == 'V':
		pylab.xlabel(xlabel + " [V]")
	elif xlabel[0] == 'I':
		pylab.xlabel(xlabel + " [A]")
	# here we hope all variables are of the same type
	if y2y1_list[0][0] == 'V':
		pylab.ylabel("V [V]")
	elif y2y1_list[0][0] == 'I':
		pylab.ylabel("I [A]")
	if log:
		ax.set_xscale('log')
		ax.set_yscale('log')
	return fig


def plot_counts_vs_wins(counts, current, todisk=None):
	counts_allwin = []
	for i in range(5):
		c = counts.get_counts_array(current, bitindex=i)
		print c
		counts_allwin.append(c)
	plotx = [numpy.arange(counts.nwin)]

	plots = zip(plotx*5, counts_allwin)

	# plot!
	fig = p.figure()
	p.hold(True)
	labels = ("Range DEF", "Range 0", "Range 1", "Range 2", "Range 3")
	fmts = ("b", "g", "r", "c", "m")
	for i in range(5):
		p.bar(plots[i][0], plots[i][1], label=labels[i], color=fmts[i])
	p.hold(False)
	p.title("Counts vs 40us windows\nIin = "+"{0:.3e}".format(current)+" A" )
	p.xlabel('#window []')
	p.ylabel('Counts []')
	ax = p.gca()
	ax.set_xlim(plots[i][0][0], plots[i][0][-1])
	p.legend()
	if todisk is not None:
		save_figure(todisk, '_vstime')
		#p.show()
	return fig

def save_figure(filename, fig):
	fig.set_size_inches(20, 10)
	pylab.savefig(filename, dpi=100, bbox_inches='tight', format=options.plotting_outtype)
	return

def plot_file(title, x, y2y1_list, filename, analysis, outfilename):
	"""Reads out data from a file and plots it.
	Deprecated, use plot_file, it passes the file directly to gnuplot.
	"""
	fig = pylab.figure()
	setup_plot(fig, title, x, y2y1_list, analysis)
	gdata = []
	gx = read_data(filename, x)
	pylab.hold(True)
	for y2label, y1label in y2y1_list:
		if y1label is not None and y1label != '':
			data1 = read_data(filename, y1label)
			ylabel = y2label+"-"+y1label
		else:
			ylabel = y2label
			data1 = 0
		data2 =  read_data(filename, y2label)
		pylab.plot(gx, data2-data1, options.plotting_style, label=ylabel+" ("+analysis+")",)
	pylab.hold(False)
	pylab.legend()

	if outfilename is not None and options.plotting_outtype is not None:
		save_figure(outfilename, fig)
	return


def show_plots():
	pylab.show()


if __name__ == '__main__':
	filename = 'colpitts_graph.tran'
	labels = read_data_header(filename)
	print labels
	for label in labels[1:]:
                print label
		plot_file("Plot "+label, labels[0], [label], [None], filename, "tran", "plot-"+label+".png")
