
# coding: utf-8

# In[3]:

from mvpa2.suite import *
from mvpa2.datasets.eventrelated import fit_event_hrf_model
import numpy as np
import os
from glob2 import glob
from os.path import join as _opj
from group_level_map import apply_warp
from datetime import datetime

import pickle


def run_searchlight(op, subjectdir, conf, output_dir,TR=2):
	mask_name = conf.mask_name
	conditions = conf.conditions_to_compare
	flavor = conf.flavor
	study_path = op.study_dir()
	subcode = subjectdir.subcode()

	for condition in conditions:
		did_run = True
		output = _opj(output_dir, '*{}{}*'.format(condition[0], condition[1]))
		if conf.num_of_permutations > 0:
			output = "{}_perm{}".format(output,conf.num_of_permutations)
		if len(glob(output)) == 0:
			did_run = False
	if did_run:
		print "already ran all sl for {}".format(output_dir)
		return

	fds = conf.get_ds(study_path, subcode, conf, mask_name, flavor, TR)
	print fds.summary()
	warp = glob(_opj(study_path,'sub{:0>3d}'.format(subcode), '**', conf.mvpa_tasks[0], 'reg', 'example_func2standard_warp.nii.gz'))[0]

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	for pair in conditions:
		permute = AttributePermutator('condition', limit='chunks')
		print conf.num_of_permutations+1
		for j in xrange(conf.num_of_permutations+1):
			prefix = conf.get_cond_prefix(pair)
			cond_ds = fds[np.array([c in pair for c in fds.sa['condition']])]
			if j > 0:
				cond_ds = permute(cond_ds)
				prefix = "{}_perm{}".format(prefix,j)
			print prefix
			output_basename = os.path.join(output_dir, prefix)
			if(len(glob(output_basename+"*")) > 0):
				print "sl already ran {}".format(j)
				continue

			kwa = {'voxel_indices': conf.get_neighbourhood_strategy(cond_ds)}
			qe = IndexQueryEngine(**kwa)
			# init the searchlight with the queryengine
			sl = Searchlight(conf.get_sl_measure(), queryengine=qe, roi_ids=None,
	                       enable_ca=['roi_sizes', 'roi_feature_ids'])
			print "starting sl {}".format(datetime.now())
			sl_map = sl(cond_ds)
			print "finished sl {}".format(datetime.now())

			pickle.dump(sl_map, open("{}_sl_map.p".format(output_basename), "wb"))
			acc_results = map2nifti(sl_map,
	                           imghdr=fds.a.imghdr)
			acc_nii_filename = '{}-acc.nii.gz'.format(output_basename)
			acc_results.to_filename(acc_nii_filename)
			#do_searchlight(cond_ds,k,os.path.join(output_dir, prefix))

			out_filename = acc_nii_filename.replace('.nii.gz', '_mni.nii.gz')
			apply_warp(acc_nii_filename, warp, out_filename)

