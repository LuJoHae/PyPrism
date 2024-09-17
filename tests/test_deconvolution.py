import unittest
import pyprism
import numpy as np

from pyprism.deconvolution import deconvolution_parallel


class TestDeconvolutionParallel(unittest.TestCase):
    def test_deconvolution_parallel(self):
        test_status = True
        for seed in range(20):
            reference = pyprism.simulation.generate_toy_single_cell_data(cell_number=4, gene_number=1000, seed=seed)
            bulk, theta = pyprism.simulation.generate_bulk_data_from_single_cell_data(reference, number_of_samples=1, scale=10, seed=seed)
            bulk, reference = pyprism.utils.reduce_to_common_var(bulk, reference)
            B = pyprism.deconvolution.deconvolution_parallel(bulk=bulk, reference=reference, n=1000)
            test_status = test_status and np.all(pyprism.deconvolution.calculate_fractions(B) - theta.X / theta.X.sum() < 0.02)
        self.assertTrue(test_status)


if __name__ == '__main__':
    unittest.main()