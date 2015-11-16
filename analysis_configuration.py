from knn_neighbourhood import KNNNeighbourhood
from mvpa2.misc.neighborhood import Sphere
from mvpa2.clfs.svm import LinearCSVMC
from mvpa2.base.node import ChainNode
from mvpa2.generators.resampling import Balancer
from mvpa2.generators.partition import NFoldPartitioner
from mvpa2.measures.base import CrossValidation
from mvpa2.misc.errorfx import mean_match_accuracy
from mvpa2.mappers.fx import mean_sample
from ds_creation import *

def get_linear_svm_measure():
		clf = LinearCSVMC(space='condition')
		splt = ChainNode([NFoldPartitioner(),Balancer(attr='condition',count=1,limit='partitions',apply_selection=True)],space='partitions')
		#splt = NFoldPartitioner()
		cv = CrossValidation(clf, splt,
								errorfx=mean_match_accuracy,
								enable_ca=['stats'], postproc=mean_sample())
		return cv
class AnalysisConfiguration(object):
	def __init__(self):
		self.study_name = 'LP'

		self.func_seg = True
		#self.standard_mask = '/home/user/data/LP/models/model002/mvpa_mask.nii.gz'
		self.flavor = 'mcf_4trim'
		self.mask_name = '_shadowreg_cluster_mask_zstat2'
		self.dir_suffix = 'linear'

		self.num_of_permutations = 50
		self.conditions_to_compare = [['G1', 'G4'], ['G2', 'G3']]

		self.mvpa_tasks = ['task001']
		self.num_of_volumes_to_delete = 4

		self.neighbourhood_type = 'knn'
		self.neighbourhood_size = 100

		self.ds_type = 'beta_trial'

	def dir_name(self):
		return '{}_{}_{}'.format(self.flavor, self.mask_name, self.dir_suffix)

	def get_cond_prefix(self,conditions):
		return "{}_{}{}_{}{}".format(self.mask_name, self.neighbourhood_type,
		                      self.neighbourhood_size, conditions[0], conditions[1])

	def get_neighbourhood_strategy(self,dataset):
		if self.neighbourhood_type == 'knn':
			return KNNNeighbourhood(self.neighbourhood_size, dataset.fa['voxel_indices'])
		elif self.neighbourhood_type == 'sphere':
			return Sphere(self.neighbourhood_size)

	def get_sl_measure(self):
		return get_linear_svm_measure()

	def get_ds(self,study_path, subcode, conf, mask_name, flavor, TR):
		if self.ds_type == 'beta_trial':
			return create_betas_per_trial_with_pymvpa(study_path, subcode, conf, mask_name, flavor, TR)
		elif self.ds_type == 'beta_run':
			return create_betas_per_run_with_pymvpa(study_path, subcode, conf, mask_name, flavor)

