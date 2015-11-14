__author__ = 'user'
from sklearn.neighbors import NearestNeighbors
import numpy as np

class KNNNeighbourhood(object):
    def __init__(self, num_neighbours, voxel_indices):
        """
        """

        self._num_neighbours = num_neighbours
        self._voxel_indices = map(tuple, voxel_indices) # we want to convert each voxel index to tuple, since we need to hash it down the road..

        nbrs = NearestNeighbors(n_neighbors=self._num_neighbours, algorithm='ball_tree').fit(voxel_indices)
	self._distances, self._indices = nbrs.kneighbors(voxel_indices)

	# for faster lookups we are creating a dictionary that will map coordinates to their index position
	# for example:
	# [(1, 1), (2, 2), (3, 3), (4, 4)] will be mapped to {(1, 1): 0, (2, 2): 1, (3, 3): 2, (4, 4): 3}
	self._coordinates_to_index = dict(zip(self._voxel_indices, range(len(self._voxel_indices))))

    def __repr__(self, prefixes=[]):
        return str(self._voxel_indices)

    def train(self, dataset):
        pass

    def __call__(self, coordinate):
        """Get all coordinates within diameter

        Parameters
        ----------
        coordinate : sequence type of length 3 with integers

        Returns
        -------
        list of tuples of size 3

        """
        # type checking
        coordinate = tuple(np.asanyarray(coordinate))
	# TODO: exception..
	coordinate_idx = self._coordinates_to_index[coordinate]
	
	neighbours = [self._voxel_indices[neighbour] for neighbour in self._indices[coordinate_idx]]

	return neighbours

