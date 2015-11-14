#!/usr/bin/python
from glob2 import glob
from os.path import join as _opj
import sys
from nipype.interfaces import fsl
import subprocess
import os
import numpy as np
import nibabel as nib

def generate_group_level_map(mni_files, output_dir):

    ext_files = filter(lambda f: '_G1G4' in f, mni_files)
    ext_prefix = "ext_{}".format(len(ext_files))
    int_files = filter(lambda f: '_G2G3' in f, mni_files)
    int_prefix = "int_{}".format(len(int_files))

    print "14:{}".format(len(ext_files))
    print "23:{}".format(len(int_files))
    thrs = [0.58, 0.6]
    ext_all, ext_avg = calc_summary_niis(ext_files, output_dir, ext_prefix)
    int_all, int_avg = calc_summary_niis(int_files, output_dir, int_prefix)
    for thr in thrs:
        ext_file = process_files(ext_prefix, output_dir, thr, ext_all, ext_avg)
        int_file = process_files(int_prefix, output_dir, thr, int_all, int_avg)
        standard_image = fsl.Info.standard_image('MNI152_T1_2mm_brain.nii.gz')
        cmd = "fslview {} {} -l Red -t 0.5 -b 3,10 {} -l Blue -t 0.5 -b 3,10 ".format(standard_image,
                                                                                       ext_file,
                                                                                       int_file)
        pro = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               shell=True, preexec_fn=os.setsid)


def process_files(prefix, output_dir, thr,all_file, avg_file):
    from scipy import ndimage
    data = avg_file.get_data()
    cluster_map, n_clusters = ndimage.label(data > thr)
    output_file = _opj(output_dir, "{}_thr_{}.nii.gz".format(prefix, thr))
    nib.save(nib.Nifti1Image(cluster_map, None, avg_file.header), output_file)
    data = all_file.get_data()
    thr_data = data > thr
    res = np.sum(thr_data, 3)
    output_file = _opj(output_dir, "{}_sum_thr_{}.nii.gz".format(prefix, thr))
    nib.save(nib.Nifti1Image(res, None, avg_file.header), output_file)

    return output_file


def calc_summary_niis(in_files, output_dir, prefix):
    all_file = _opj(output_dir, '{}_all.nii.gz'.format(prefix))
    avg_file = _opj(output_dir, '{}_avg.nii.gz'.format(prefix))
    merge = fsl.Merge(in_files=in_files,
                      dimension='t',
                      merged_file=all_file)
    merge.run()
    mean = fsl.maths.MeanImage(in_file=all_file, dimension='T', out_file=avg_file)
    mean.run()
    all_nii = nib.load(all_file)
    avg_nii = nib.load(avg_file)
    all_data = all_nii.get_data()
    med_data = np.median(all_data,3)
    output_file = _opj(output_dir, "{}_median.nii.gz".format(prefix))
    nib.save(nib.Nifti1Image(med_data, None, avg_nii.header), output_file)
    return all_nii, avg_nii


def apply_warp(in_file, warp_file, out_file):
    standard_image = fsl.Info.standard_image('MNI152_T1_2mm_brain.nii.gz')
    apply_warp = fsl.preprocess.ApplyWarp(ref_file=standard_image,
                                          in_file=in_file,
                                          field_file=warp_file,
                                          interp='nn',
                                          out_file=out_file)
    apply_warp.run()


if __name__ == '__main__':
    path = sys.argv[1]
    classifier = sys.argv[2]
    output_dir = sys.argv[3]
    files = sorted(glob(_opj(path, '**', 'sub*{}*acc*.nii.gz'.format(classifier))))
    generate_group_level_map(path, files, classifier, output_dir)
