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

def setup_plot(fig, title, xvu, yvu, log=False, xlog=False, ylog=False):
	""" fig: the figure
	title: plot title:
	xvu: tuple defined as xvu = (xvar, xunit)
	yvu: list of tuples defined as yvu += [(yvarN, yunitN)]

	returns: fig
	"""
	#xvar, xunit = xvu
	pylab.title(title.upper())
	ax = pylab.gca()
	pylab.xlabel('%s [%s]' % xvu)
	yunits = []
	yinitials = []
	for yv, yu in yvu:
		if not yu in yunits:
			yunits.append(yu)
			yinitials.append(yv[0])
	ylabel = ""
	for yi, yu in zip(yinitials, yunits):
		ylabel += "%s [%s] / " % (yi, yu)
	ylabel = ylabel[:-3]
	pylab.ylabel(ylabel)
	if log or xlog:
		ax.set_xscale('log')
	if log or ylog:
		ax.set_yscale('log')
	return fig

def save_figure(filename, fig):
	fig.set_size_inches(20, 10)
	pylab.savefig(filename, dpi=100, bbox_inches='tight', format=options.plotting_outtype)
	return

def plot_results(title, xvarname, y2y1_list, results, outfilename):
	"""Plot the results.
	"""
	fig = pylab.figure()
	analysis = results.get_type().upper()
	xunit = results.units[xvarname]
	gdata = []
	x = results[xvarname]
	yvu = []

	for y2label, y1label in y2y1_list:
		if y1label is not None and y1label != '':
			data1 = results[y1label]
			line_label = y2label+"-"+y1label
		else:
			line_label = y2label
			data1 = 0
		data2 =  results[y2label]
		yvu += [(line_label, results.units[y2label])]
		gdata.append((data2-data1, line_label))

	if xvarname == 'w':
		xlog = True
	else:
		xlog = False
	setup_plot(fig, title, (xvarname, xunit), yvu, xlog=xlog)

	pylab.hold(True)
	for y, label in gdata:
		# pylab wants matrices in the form: N,1, while results come as (1, N) -> Transpose
		pylab.plot(x.T, y.T, options.plotting_style, label=label+" ("+analysis+")",)
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
