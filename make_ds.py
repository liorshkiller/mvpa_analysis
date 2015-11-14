#!/usr/bin/python

import sys
from os.path import join as _opj
import os
from mvpa2.datasets.sources.openfmri import OpenFMRIDataset
from mvpa2.datasets.eventrelated import fit_event_hrf_model
from mvpa2.base.hdf5 import h5save
from mvpa2.mappers.detrend import poly_detrend
from mvpa2.mappers.zscore import zscore
# from nilearn.image import smooth_img
import nibabel as nb
import numpy as np
from SubjectDir import SubjectDir
from OpenFMRIData import OpenFMRIData
from mvpa2.base.dataset import vstack
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
    poly_detrend(ds, polyord=2, chunks_attr='chunks')
    zscore(ds, chunks_attr='chunks', dtype='float32')
    return ds


def make_ds(subject_dir, datapath, flavor, mask_file, dir_name=None ):
    if not dir_name:
        dir_name = flavor
    result_dir = _opj(subject_dir.model_dir(1), dir_name, 'ds')
    if not os.path.isdir(result_dir):
        os.makedirs(result_dir)
    sub = subject_dir.subcode()
    ds14_filename = _opj(result_dir, 'sub%.3i_14_hrf.hdf5' % sub)
    ds23_filename = _opj(result_dir, 'sub%.3i_23_hrf.hdf5' % sub)
    ds_int_ext_filename = _opj(result_dir, 'sub%.3i_int_ext_hrf.hdf5' % sub)
    if os.path.isfile(ds_int_ext_filename):
        print "ds exists"
        return [ds14_filename, ds23_filename]
    of = OpenFMRIDataset(datapath)

    ds = of.get_model_bold_dataset(
    model_id=1, subj_id=sub,
    flavor=flavor,
    mask= _opj(
        datapath, 'sub%.3i' % sub, 'masks', 'task001',mask_file), #'grey.nii.gz'),
    # preproc_img=smooth,
    preproc_ds=detrend,
    modelfx=fit_event_hrf_model,
    time_attr='time_coords',
    condition_attr='condition')

    ds14 = ds[np.array([c in ['G1', 'G4'] for c in ds.sa['condition']])]
    ds23 = ds[np.array([c in ['G2', 'G3'] for c in ds.sa['condition']])]

    print "{:0>3d}-ds14 {},{}".format(sub, ds14.shape, ds14.sa.condition)
    print "{:0>3d}-ds23 {},{}".format(sub, ds23.shape, ds23.sa.condition)
    h5save(ds14_filename, ds14)
    h5save(ds23_filename, ds23)
    ds14.sa.condition = ['ext'] * len(ds14.sa.condition)
    ds23.sa.condition = ['int'] * len(ds23.sa.condition)
    int_ext_ds = vstack([ds14, ds23],0)
    int_ext_ds.a.imghdr = ds.a.imghdr
    h5save(ds_int_ext_filename, int_ext_ds)
    print "{:0>3d}-ds23 {},{}".format(sub, int_ext_ds.shape, int_ext_ds.sa.condition)
    return [ds14_filename, ds23_filename]

def main():
    sub = int(sys.argv[1])
    data_dir = os.environ.get('DATA_DIR') or '/home/user/data'
    study_name = os.environ.get('STUDY_NAME') or 'LP'
    flavor = 'mcf'
    op = OpenFMRIData(data_dir, study_name)
    make_ds(op.subject_dir(subcode=sub), op.study_dir(), flavor)


if __name__ == '__main__':
    main()
