import os
import sys
import getopt
import numpy as np
import scipy
import scipy.io
import re

try:
	opts, args = getopt.getopt(sys.argv[1:], '', ['alg=', 'max_cluster=', 'dataset=', 'width=', 'height=', 'ldis='])
except getopt.GetoptError:
	print 'wrong parameters'
	sys.exit(1)

opts = dict(opts)
if len(args) == 1:
	input_file = args[0]
elif len(args) == 0:
	print 'No input file given'
	sys.exit(1)
else:
	print 'More than 1 input files are given'
	sys.exit(1)
if '--alg' in opts:
	alg = int(opts['--alg'])
else:
	print 'Please specify an algorithm in the mat-file'
	sys.exit(1)
if '--max_cluster' in opts:
	max_cluster = int(opts['--max_cluster'])
else:
	max_cluster = -1
if '--dataset' in opts:
	dataset_name = opts['--dataset']
else:
	dataset_name = 'Data set'
if '--width' in opts:
	width = opts['--width']
else:
	width = '11in'
if '--height' in opts:
	height = opts['--height']
else:
	height = '8.5in'
if '--ldis' in opts:
	ldis = opts['--ldis']
else:
	ldis = '6cm'

v_reg = '\^'
v_reg = re.compile(v_reg)
v2_reg = '#'
v2_reg = re.compile(v2_reg)
v3_reg = '&'
v3_reg = re.compile(v3_reg)
exp_reg = '_'
exp_reg = re.compile(exp_reg)

class Node:
	def __init__(self, isroot, count, voc=None, percent=None, priority=None):
		self.isroot = isroot
		self.count = count
		self.voc = voc
		self.percent = percent
		self.priority = priority

	def latex(self, f, indent):
		indent += 2
		if self.isroot:
			f.write(' '*indent + r'\Tree [.{\Large{\textbf{' + dataset_name + ', Total: ' + str(self.count) + '}}}\n')
		else:
			if hasattr(self, 'left_child'):
				color = 'blue'
			else:
				color = 'red'
			f.write(' '*indent + r'.\node[GenericNodeStyle, draw=' + color + '] {\n')
			indent += 2
			f.write(' '*indent + r'\begin{tabular}{c}' + '\n')
			f.write(' '*indent + r'Count: ' + str(self.count) + r' \\' + '\n')
			line_num = 0
			while line_num < 10 and (line_num < len(self.voc) and self.voc[line_num].shape != (1, 0)):
				v = unicode.encode(self.voc[line_num][0], 'utf-8')
				v = v_reg.sub('\\^', v)
				v = v2_reg.sub('\\#', v)
				v = v3_reg.sub('\\&', v)
				v = exp_reg.sub('\\_', v)
				f.write(' '*indent + "`" + v + "'" + r' \\' + '\n')
				line_num += 1
			f.write(' '*indent + r'\end{tabular}' + '\n')
			indent -= 2
			f.write(' '*indent + '};\n')
		if hasattr(self, 'left_child'):
			indent += 2
			f.write(' '*indent + '[\n')
			self.left_child.latex(f, indent)
			f.write(' '*indent + ']\n')
			f.write(' '*indent + '[\n')
			self.right_child.latex(f, indent)
			f.write(' '*indent + ']\n')
			indent -= 2
		if self.isroot:
			f.write(' '*indent + ']\n')
		indent -= 2

def insert(root, old_node, new_node, new_label):
	split_node = root
	mask = 1
	level = 0
	while hasattr(split_node, 'left_child'):
		if mask & new_label:
			split_node = split_node.right_child
		else:
			split_node = split_node.left_child
		mask <<= 1
		level += 1
	split_node.left_child = old_node
	split_node.right_child = new_node
	return level

def get_old_label(new_label):
	mask = 1
	while mask <= new_label:
		mask <<= 1
	return new_label - (mask>>1)

def map_label(unique_labels):
	label_map = {}
	idx = 0
	for i in unique_labels:
		label_map[i] = idx
		idx += 1
	return label_map

data = scipy.io.matlab.mio.loadmat(input_file)
name, ext = os.path.splitext(input_file)

clusters = data['final_cluster_summary'][alg-1]
counts = data['final_count_summary'][alg-1]
vocs = data['final_voc_summary'][alg-1]
percents = data['final_percent_summary'][alg-1]
priorities = data['final_priority_summary'][alg-1]
if max_cluster > 0:
	max_cluster -= 1
	clusters = clusters[0:max_cluster]
	counts = counts[0:max_cluster]
	vocs = vocs[0:max_cluster]
	percents = percents[0:max_cluster]
	priorities = priorities[0:max_cluster]
num_splits = len(clusters)

cur_labels = np.array([0])
max_level = 0
for one_split in range(num_splits):
	cluster = clusters[one_split]
	count = counts[one_split]
	voc = vocs[one_split]
	percent = percents[one_split]
	priority = priorities[one_split]
	unique_labels = np.unique(cluster)
	unique_labels = np.setdiff1d(unique_labels, np.array([-1]))
	new_label = np.setdiff1d(unique_labels, cur_labels)[0]
	label_map = map_label(unique_labels)
	old_label = get_old_label(new_label)
	old_idx = label_map[old_label]
	new_idx = label_map[new_label]
	if one_split == 0:
		root = Node(True, int(count[old_idx][0]) + int(count[new_idx][0]))
	old_node = Node(False, count[old_idx][0], voc[:,old_idx], None, priority[0][old_idx])
	new_node = Node(False, count[new_idx][0], voc[:,new_idx], None, priority[0][new_idx])
	level = insert(root, old_node, new_node, new_label)
	if level > max_level:
		max_level = level
	cur_labels = unique_labels

f_output = open(name + '.tex', 'w')
f_output.write('\\documentclass{article}\n')
f_output.write('\\usepackage[paperwidth=' + width + ', paperheight=' + height + ', left=1in, right=1in, top=1in, bottom=1in]{geometry}\n')
f_output.write('''
\\usepackage{tikz}
\\usepackage{tikz-qtree}
\\usepackage{multirow}

\\tikzset
{
  GenericNodeStyle/.style =
  {
    shape = rectangle,
    rounded corners = 2mm,
    scale = 1.0,
    thick,
  }
}

\\begin{document}
\\thispagestyle{empty}

\\begin{center}
''')
f_output.write('\\begin{tikzpicture}[level distance=' + ldis + ']\n')

root.latex(f_output, 0)

f_output.write(''';
\\end{tikzpicture}
\\end{center}

\\end{document}
''')
f_output.close()

os.system('pdflatex ' + name + '.tex')
os.system('evince ' + name + '.pdf')
