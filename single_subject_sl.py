#!/usr/bin/python

from mvpa2.measures.searchlight import Searchlight
from mvpa2.clfs.svm import LinearCSVMC
from mvpa2.generators.partition import NFoldPartitioner
from mvpa2.measures.base import CrossValidation
from mvpa2.misc.neighborhood import IndexQueryEngine
from mvpa2.misc.errorfx import mean_match_accuracy
from mvpa2.mappers.fx import mean_sample
from mvpa2.cmdline.helpers import arg2ds
from mvpa2.generators.permutation import AttributePermutator
from mvpa2.clfs.stats import MCNullDist
from mvpa2.datasets.mri import map2nifti
from knn_neighbourhood import KNNNeighbourhood
from mvpa2.base.node import ChainNode
from mvpa2.generators.resampling import Balancer
import sys
import os.path
from glob2 import glob
"""
From: https://github.com/PyMVPA/PyMVPA/blob/master/mvpa2/cmdline/cmd_searchlight.py
"""


def _fill_in_scattered_results(sl, dataset, roi_ids, results):
    """this requires the searchlight conditional attribute 'roi_feature_ids'
    to be enabled"""
    import numpy as np
    from mvpa2.datasets import Dataset

    resmap = None
    probmap = None
    for resblock in results:
        for res in resblock:
            if resmap is None:
                # prepare the result container
                resmap = np.zeros((len(res), dataset.nfeatures),
                                  dtype=res.samples.dtype)
                if 'null_prob' in res.fa:
                    # initialize the prob map also with zeroes, as p=0 can never
                    # happen as an empirical result
                    probmap = np.zeros((dataset.nfeatures,) + res.fa.null_prob.shape[1:],
                                       dtype=res.samples.dtype)
                observ_counter = np.zeros(dataset.nfeatures, dtype=int)
            # project the result onto all features -- love broadcasting!
            #print "averaging"
            resmap[:, res.a.roi_feature_ids] += res.samples
            if not probmap is None:
                probmap[res.a.roi_feature_ids] += res.fa.null_prob
            # increment observation counter for all relevant features
            observ_counter[res.a.roi_feature_ids] += 1
    # when all results have been added up average them according to the number
    # of observations
    observ_mask = observ_counter > 0
    resmap[:, observ_mask] /= observ_counter[observ_mask]
    result_ds = Dataset(resmap,
                        fa={'observations': observ_counter})
    if not probmap is None:
        # transpose to make broadcasting work -- creates a view, so in-place
        # modification still does the job
        probmap.T[:, observ_mask] /= observ_counter[observ_mask]
        result_ds.fa['null_prob'] = probmap.squeeze()
    if 'mapper' in dataset.a:
        import copy
        result_ds.a['mapper'] = copy.copy(dataset.a.mapper)
    return result_ds


def do_searchlight(glm_dataset, radius, output_basename, with_null_prob=False, clf=LinearCSVMC(space='condition')):
    if(len(glob(output_basename+"*")) > 0):
        print "sl already ran"
        return
    splt = ChainNode([NFoldPartitioner(),Balancer(attr='condition',count=1,limit='partitions',apply_selection=True)],space='partitions')
    #splt = NFoldPartitioner()
    cv = CrossValidation(clf, splt,
                         errorfx=mean_match_accuracy,
                         enable_ca=['stats'], postproc=mean_sample())
    distr_est = []
    if with_null_prob:
        permutator = AttributePermutator('condition', count=100,
                                         limit='chunks')
        distr_est = MCNullDist(permutator, tail='left',
                               enable_ca=['dist_samples'])
        """
        repeater   = Repeater(count=100)
        permutator = AttributePermutator('condition', limit={'partitions': 1}, count=1)
        null_cv = CrossValidation(clf, ChainNode([splt, permutator],space=splt.get_space()),
                      postproc=mean_sample())
        null_sl = sphere_searchlight(null_cv, radius=radius, space='voxel_indices',
                         enable_ca=['roi_sizes'])
        distr_est = MCNullDist(repeater,tail='left', measure=null_sl,
                       enable_ca=['dist_samples'])

        sl = sphere_searchlight(cv, radius=radius, space='voxel_indices',
                                null_dist=distr_est,
                                enable_ca=['roi_sizes', 'roi_feature_ids']
                                # ,result_fx = _fill_in_scattered_results # average across all spheres
                                )
        """
    else:
        kwa = {'voxel_indices': KNNNeighbourhood(radius, glm_dataset.fa['voxel_indices'])}
        qe = IndexQueryEngine(**kwa)
        # init the searchlight with the queryengine
        sl = Searchlight(cv, queryengine=qe, roi_ids=None,
                       enable_ca=['roi_sizes', 'roi_feature_ids']
                       # ,results_fx = _fill_in_scattered_results # average across all spheres
                          )
       #;v sl = sphere_searchlight(cv, radius=radius, space='voxel_indices',

                                # ,result_fx = _fill_in_scattered_results # average across all spheres
                               # )
    # ds = glm_dataset.copy(deep=False,
    #		       sa=['condition','chunks'],
    #		       fa=['voxel_indices'],
    #		       a=['mapper'])
    from datetime import datetime
    print "starting sl {}".format(datetime.now())
    sl_map = sl(glm_dataset)
    print "finished sl {}".format(datetime.now())
    import pickle
    pickle.dump(sl_map, open("{}_sl_map.p".format(output_basename), "wb"))
#    pickle.dump(sl.ca.roi_feature_ids, open("{}_sl_feature_ids.p".format(output_basename), "wb"))
#    print len(sl.ca.roi_feature_ids[0])
    acc_results = map2nifti(sl_map,
                           imghdr=glm_dataset.a.imghdr)
    acc_nii_filename = '{}-acc.nii.gz'.format(output_basename)
    acc_results.to_filename(acc_nii_filename)
    sl_map.samples *= -1
    sl_map.samples += 1
    niftiresults = map2nifti(sl_map,
                             imghdr=glm_dataset.a.imghdr)
    niftiresults.to_filename('{}-err.nii.gz'.format(output_basename))
    # TODO: check p value map
    if with_null_prob:
        nullt_results = map2nifti(sl_map, data=sl.ca.null_t,
                                  imghdr=glm_dataset.a.imghdr)
        nullt_results.to_filename('{}-t.nii.gz'.format(output_basename))
        nullprob_results = map2nifti(sl_map, data=sl.ca.null_prob,
                                     imghdr=glm_dataset.a.imghdr)
        nullprob_results.to_filename('{}-prob.nii.gz'.format(output_basename))
        nullprob_results = map2nifti(sl_map, data=distr_est.cdf(sl_map.samples),
                                     imghdr=glm_dataset.a.imghdr)
        nullprob_results.to_filename('{}-cdf.nii.gz'.format(output_basename))
    return sl_map


if __name__ == '__main__':
    filename = sys.argv[1]
    radius = int(sys.argv[2])
    print filename
    output_basename = os.path.join('{}_r{}_c-{}'.format(filename, radius, 'linear'))
    print output_basename
    ds = arg2ds([filename])
    do_searchlight(ds, radius, output_basename,
                   False)
