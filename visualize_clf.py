#!/usr/bin/python
from matplotlib import pyplot as plt
import numpy as np
import pickle
from glob2 import glob
from scipy import ndimage
import sys


def visualize_clf(file_path):
	ext_pattern = "14"
	int_pattern = "23"
	path = "{}/**/*{}*.p".format(file_path,ext_pattern)
	
	files = glob(path)
	print files
	thresholds = np.arange(0.65,1,0.05)
	file_dict = dict()
	for f in files:
		filename = f[f.rfind('/')+1:]
		sub = filename[:filename.find('_')]
		pair = (f,f.replace(ext_pattern,int_pattern))
		print pair
		if sub in file_dict:
			file_dict[sub].append(pair) 
		else:
			file_dict[sub]=[pair]
	print file_dict
	for sub,file_list in file_dict.iteritems():
		fig = plt.figure()	
		cell_text = []
		col_labels= []
		file_list = sorted(file_list)
		for i,pair in enumerate(file_list):
			print pair
			f = pair[0]
			sl = pickle.load(open(f,'rb'))
			data = sl.samples[0]
			fig.add_subplot(4,4,i+1)
			title = f[f.find('-')+1:]
			plt.title(title)
			col_labels.append(title)
			plt.hist(data)
			coltext = []
			print title
			for thr in thresholds:
			    data_3d = sl.a.mapper.reverse1(sl.samples)[0]
			    cluster_map, n_clusters = ndimage.label(data_3d > thr)
			    cluster_sizes = np.bincount(cluster_map.ravel())[1:]
			    if len(cluster_sizes) != 0:
			        coltext.append("{}".format(np.max(cluster_sizes)))
			    else:
				coltext.append(0)
			cell_text.append(coltext)
		ax = fig.add_subplot(4,4,len(files)+2)
		ax.axis('off')
		print len(cell_text)
		plt.table(cellText= cell_text,rowLabels=col_labels, 
				colLabels=thresholds,loc='center right')
		plt.savefig('{}.png'.format(sub))
if __name__ == '__main__':
	path = str(sys.argv[1])
	visualize_clf(path)
	raw_input()

