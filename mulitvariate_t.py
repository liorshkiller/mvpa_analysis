import numpy as np
from mvpa2.measures.base import Measure
from mvpa2.datasets.base import Dataset
from mvpa2.base.param import Parameter
import numpy.matlib as npm

class MultiT(Measure):
	"""
	"""

	is_trained = True # Indicate that this measure is always trained.
	# example of a parameter
	#    center_data = Parameter(False, constraints='bool', doc="""\
	#          If True then center each column of the data matrix by subtracing the
	#          column mean from each element. This is recommended especially when
	#          using pairwise_metric='correlation'.""")

	def __init__(self, **kwargs):
		Measure.__init__(self, **kwargs)

	def printMatrix(mat,text):
		print(text + "\n")
		print(mat)
		print("\n")

	def _checkDeltaForExceptions(self,delta):
		if np.size(delta,0) <= 3:
			raise NameError('You have less than 4 trials. Multi-T needs a min. of 3 trials / observations')
		if (~delta.any(axis=0)).any():
			raise NameError('You have one or more column (feature) that has all zeros ')

	def _checkLabelsForExceptions(self,labels):
		if np.size(labels,0)<=6:
			raise NameError('You do not have enough trials (observations), the min. is 3 from each class')
		# check that you have a balanced classes, and at least 3 observations / class
		uniqlabels = np.unique(labels)
		if np.size(uniqlabels) > 2:
			raise NameError('Currently we only support 2 classes')
		# check that oyu have equal number of labels
		countlabelsA = np.size(np.where(uniqlabels[0] == labels),0)
		countlabelsB = np.size(np.where(uniqlabels[1] == labels),0)

		if countlabelsA != countlabelsB:
			raise NameError('You do not have an equal number of labels from each class')


	def calcmultit(self,data, labels):
		"""
		This function recieves as input an array and calculates the multivariate-T value as described in:
		Srivastava, M. S. and M. Du (2008). "A test for the mean vector with fewer observations than the dimension."
		Journal of Multivariate Analysis 99(3): 386-402.
		@param array data: an array of numbers, with at least 6 rows (3 trials), no zero columns
		@param array labels: an array of trial labels, must be balanced (equal number of trials labels from class A and B)
		:return:
		"""
		# initial check to make sure that data and labels fit
		if np.size(labels,0) != np.size(data,0):
			raise NameError('The length of labels isn''t equal to the number of rows in data, check your inputs')

		#check labels for exception
		self._checkLabelsForExceptions(labels)

		#get unique labels
		uniqlabels = np.unique(labels)
		idxlabelsA = np.where(uniqlabels[0] == labels)[0]
		idxlabelsB = np.where(uniqlabels[1] == labels)[0]

		# ccompute delta of data
		delta = data[idxlabelsA,:] - data[idxlabelsB,:]
		if (~delta.any(axis=0)).any():
			return 0
		# check delta for exceptions:
		self._checkDeltaForExceptions(delta)

		# calc N,n and p
		N = np.size(delta,0)
		# printMatrix(N,"N is:")
		n = N-1
		# printMatrix(n,"n is:")
		p = np.size(delta,1)
		# printMatrix(p,"p is:")


		# calc cov matrix (transpose to get it like matlab output)
		covDelta = np.cov(delta.T)
		# printMatrix(covDelta,'cov delta:')

		# calc identity matrix of cov delta
		eyeCovDelta = np.eye(np.size(covDelta,0),np.size(covDelta,1))
		# printMatrix(eyeCovDelta,"eye cov delta:\n")

		# 1/ sqrt(diag) of cov delta
		oneDivSqrtDiagCovDelta = 1/(np.sqrt(np.diag(covDelta)))
		# printMatrix(oneDivSqrtDiagCovDelta,"oneDivSqrtDiagCovDelta \n")

		# rep mat of cov delta
		repMatOfOneDivSqrtDiagCovDelta = npm.repmat(oneDivSqrtDiagCovDelta,np.size(covDelta,1),1).T
		# printMatrix(repMatOfOneDivSqrtDiagCovDelta,"rep mat of one div sqrt cov data")

		# calc cov delta sqrt
		covDeltaSqrt = eyeCovDelta * repMatOfOneDivSqrtDiagCovDelta
		# printMatrix(covDeltaSqrt,"cov detla sqrt")

		# calc mean of delta
		meanDelta = np.average(delta,0)

		# cal corr
		corrDelta = (covDeltaSqrt.dot(covDelta)).dot(covDeltaSqrt)

		# calc trace r2 of cor delta
		traceR2 = np.trace((corrDelta.T).dot(corrDelta))
		# printMatrix(traceR2,'trace r2')

		# calc cov delta diag thingi
		covDeltaDiag = (1/covDelta) * eyeCovDelta


		####
		## calc multi t
		####
		# numerator
		numeratorT = ((N * meanDelta).dot(covDeltaDiag)).dot(meanDelta.T) - ((n*p)/(n-2))
		# print(numeratorT)
		# denominator
		denmoniator = 2 * (traceR2 - (float(p**2)/float(n)) )
		# print(denmoniator)
		# cpn fix
		cPn = 1 + (traceR2 / (pow(p,(float(3)/float(2)))) )
		# print(cPn)
		multiT = (float(numeratorT)) / np.sqrt(float(denmoniator) *float(cPn))
		# print(multiT)
		return multiT

	def _call(self,ds):
		if 'id' in ds.sa.keys():
			dt = np.dtype([('data', np.float64, (ds.samples.shape[1],)), ('labels', 'S10'), ('ids', 'S10')])
			array = np.array(zip(ds.samples, ds.sa.condition, ds.sa.id),dtype=dt)
			array = np.sort(array, order='ids')
			data = array['data']
			labels = array['labels']
		else:
			data = ds.samples
			labels = ds.sa.condition
		# center data if specified
		#if self.params.center_data:
		#	data = data - np.mean(data,0)

		t = self.calcmultit(data, labels)

		# if square return value make dsm square
			# re-ad d the sample attributes -- should still be valid
		out = Dataset([t],  fa={'metrics': ['t']})

		return out