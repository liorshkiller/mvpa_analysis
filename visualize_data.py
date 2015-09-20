#!/usr/bin/python

import sys
from os.path import join as _opj
import os
from mvpa2.datasets.sources.openfmri import OpenFMRIDataset
from mvpa2.datasets.eventrelated import fit_event_hrf_model
from mvpa2.base.hdf5 import h5save
from mvpa2.mappers.detrend import poly_detrend
from mvpa2.mappers.zscore import zscore
from mvpa2.misc.plot.base import plot_samples_distance
#from nilearn.image import smooth_img
import nibabel as nb
import numpy as np
import matplotlib.pyplot as pl
"""
def smooth(img):
    # we need to preserve the original header because the smoothing function
    # fucks the TR up
    nimg = smooth_img(img, fwhm=2.0)
    return nb.Nifti1Image(nimg.get_data(),
                          img.get_affine(),
                          header=img.get_header())
"""
def detrend(ds):
	#print ds.summary()
	ds.samples = ds.samples.astype('float')
	pl.figure()
	pl.subplot(221)
	plot_samples_distance(ds, sortbyattr='chunks')
	#plot_samples_distance(ds)
	pl.title('Sample distances (sorted by chunks)')
	poly_detrend(ds, polyord=2, chunks_attr='chunks')
	pl.subplot(222)
	plot_samples_distance(ds, sortbyattr='chunks')
	pl.show()
	zscore(ds, chunks_attr='chunks', dtype='float32')
	pl.subplot(223)
	plot_samples_distance(ds, sortbyattr='chunks')
	pl.subplot(224)
#	plot_samples_distance(ds, sortbyattr='targets')
	pl.title('Sample distances (sorted by condition)')
	pl.show()
	#poly_detrend(ds, polyord=1, chunks_attr='chunks')
	#zscore(ds, chunks_attr='chunks', dtype='float32')
	return ds
def make_ds(sub, datapath, flavor):
	of = OpenFMRIDataset(datapath)
	ds = of.get_model_bold_dataset(
	    model_id=1, subj_id=sub,
	#ds = of.get_bold_run_dataset(1,1,1,
	    
	    flavor=flavor,
	    mask=_opj(
		datapath, 'sub%.3i' % sub, 'masks', 'task001_run001',
		'grey.nii.gz'),
	    #preproc_img=smooth,
	    #preproc_ds = detrend, 
	    #modelfx=fit_event_hrf_model,
	    time_attr='time_coords',
	    condition_attr='condition')
	for i in np.unique(ds.chunks)[5:]:
		detrend(ds[ds.chunks == i])

def main():
	sub = int(sys.argv[1]) 
	data_dir   	= os.environ.get('DATA_DIR') or '/home/user/data'
	study_name 	= os.environ.get('STUDY_NAME') or 'Pilot2'
	flavor		= 'mcf'
	make_ds(sub, _opj(data_dir,study_name),flavor)


if __name__ == '__main__':
	main()
	raw_input()
