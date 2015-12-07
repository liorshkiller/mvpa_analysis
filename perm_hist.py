#!/usr/bin/python
import nibabel as nib
from matplotlib import pyplot as plt
import pickle
import sys
from glob2 import glob
from os.path import join as _opj
from analysis_configuration import AnalysisConfiguration
import numpy as np
#from mvpa2.suite import *
from mvpa2.base.dataset import vstack
from mvpa2.datasets.mri import map2nifti
from scipy.stats import percentileofscore
import os
from nipype.interfaces import fsl


def apply_warp(sub_dir, in_file, out_file):
#    warp_file = _opj(sub_dir,'anatomy','reg','highres2standard_warp.nii.gz')
#    pre_mat_file = _opj(sub_dir,'BOLD','task001','reg','example_func2highres.mat')
    warp_file = _opj(sub_dir,'BOLD','task001','reg','example_func2standard_warp.nii.gz')
    standard_image = fsl.Info.standard_image('MNI152_T1_2mm_brain.nii.gz')
    apply_warp = fsl.preprocess.ApplyWarp(ref_file=standard_image,
                                          in_file=in_file,
                                          field_file=warp_file,
                                          #premat=pre_mat_file,
                                          interp='trilinear',
                                          out_file=out_file)
    apply_warp.run()
	# /usr/share/fsl/5.0/bin/applywarp --ref=reg/standard --in=stats/cope2 --out=frgrot_Qb7uSq --warp=reg/highres2standard_warp --premat=reg/example_func2highres.mat --interp=trilinear


def perm_hist(subj):
	conf = AnalysisConfiguration()
	data_dir = os.environ.get('DATA_DIR') or '/home/user/data'
	sub_dir = _opj(data_dir,conf.study_name,'sub{:0>3d}'.format(subj))
	directory = _opj(data_dir,'LP/sub{:0>3d}/results/'.format(subj))
	for pair in conf.conditions_to_compare:
			files = sorted(glob(_opj(directory,conf.dir_name(),'{}*{}{}*.p'.format(conf.mask_name,pair[0],pair[1]))))
			plt.figure()
			plt.subplot(211)
			plt.title('sub{:0>3d}-{}{}'.format(subj,pair[0],pair[1]))
			print pair, " ", len(files)
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
			perms = vstack(all_maps)
			real_f = files[-1]
			f_h = file(real_f,'r')
			real_map = pickle.load(f_h)
			color = 'crimson'
			line_width = 2
			plt.hist(np.transpose(real_map),bins=20,histtype='step',color=[color], lw = line_width)
			percentiles = np.zeros((1,len(real_map.samples[0])))
			for i,vox in enumerate(real_map.samples[0]):
			    percentiles[0,i]=percentileofscore(perms[:,i].samples.flat,vox)
			plt.subplot(212)
			print len(percentiles[0])
			plt.hist(percentiles[0],bins=20,histtype='step')
			real_map.samples=percentiles
			nii = real_f.replace("_sl_map.p", "-acc.nii.gz")
			nii_file = nib.load(nii)
			perc_results = map2nifti(real_map, imghdr=nii_file.header)
			perc_nii_filename =real_f.replace("_sl_map.p", "-percentiles_sub{:0>3d}.nii.gz".format(subj))
			perc_results.to_filename(perc_nii_filename)
			thr_prc_filename = perc_nii_filename.replace(".nii.gz","_p0.01.nii.gz")
			thr = fsl.maths.Threshold(in_file=perc_nii_filename, thresh=100,
						  out_file=thr_prc_filename)
			thr.run()
			mni_thr_filename = thr_prc_filename.replace(".nii.gz","_mni.nii.gz")
			apply_warp(sub_dir,thr_prc_filename, mni_thr_filename)

			
	plt.show()
	#plt.savefig('/tmp/sub{:0>3d}_{}{}'.format(subj,pair[0],pair[1]))
	raw_input()


if __name__ == '__main__':
	perm_hist(2)
#	for i in xrange(1,8):
#		print "sub",i
#		perm_hist(i)
