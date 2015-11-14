#!/usr/bin/python
import nibabel as nib
from matplotlib import pyplot as plt
import pickle
import sys
from glob2 import glob
from os.path import join as _opj
from analysis_configuration import AnalysisConfiguration
import numpy as np
from mvpa2.suite import *
from scipy.stats import percentileofscore

if __name__ == '__main__':
	
	conf = AnalysisConfiguration()
	subj = int(sys.argv[1])
	directory = '/home/user/data/LP/sub{:0>3d}/results/'.format(subj)
	for pair in conf.conditions_to_compare:
			files = sorted(glob(_opj(directory,'**','{}*{}{}*.p'.format(conf.mask_name,pair[0],pair[1]))))
			plt.figure()
			plt.subplot(211)
			plt.title('sub{:0>3d}-{}{}'.format(subj,pair[0],pair[1]))
			print len(files)
			all_maps = []
			for f in files[:-1]:
				f_h = file(f,'r')
				m = pickle.load(f_h)
				all_maps.append(m)
				if 'perm' in f:
					color = 'black'
					line_width = 1
				else:
					color = 'crimson'
					line_width = 2
				plt.hist(np.transpose(m),bins=20,histtype='step',color=[color], lw = line_width)
				plt.show()
			perms = vstack(all_maps)
			real_f = files[-1]
			f_h = file(real_f,'r')
			real_map = pickle.load(f_h)
			percentiles = np.zeros((1,len(real_map.samples[0])))
			for i,vox in enumerate(real_map.samples[0]):
			    percentiles[0,i]=percentileofscore(perms[:,i].samples.flat,vox)
			plt.subplot(212)
			print percentiles[0]
			plt.hist(percentiles[0],bins=20,histtype='step')
			real_map.samples=percentiles
			nii = real_f.replace("_sl_map.p", "-acc.nii.gz")
			nii_file = nib.load(nii)
			perc_results = map2nifti(real_map, imghdr=nii_file.header)
			perc_nii_filename = real_f.replace("_sl_map.p", "-percentiles.nii.gz")
			perc_results.to_filename(perc_nii_filename)
	plt.savefig('/tmp/sub{:0>3d}_{}{}'.format(subj,pair[0],pair[1]))
	raw_input()
