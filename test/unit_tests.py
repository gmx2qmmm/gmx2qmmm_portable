""" 
Run me with 
python -m unittest -v unit_tests.py
"""
import unittest
import numpy as np
from gmx2qmmm.generators.geometry import make_xyzq

class TestGeometryUtils(unittest.TestCase):

    def test_make_xyzq(self):
        
        # input data
        coords  = [0.123, 0.456, 0.789, 0.321, 0.654, 0.987]
        charges = [-0.57, +1.57]
        ref_xyzq = np.array([[0.123, 0.456, 0.789, -0.57],[0.321, 0.654, 0.987, +1.57]])

        # run the function
        res_xyzq = make_xyzq(coords,charges)

        # check it
        are_equal = np.allclose(res_xyzq, ref_xyzq)
        self.assertTrue(are_equal)




class TestProductionLevel(unittest.TestCase):
    # Production / Integration tests to be implemented by Alina
    pass



if __name__ == '__main__':
    unittest.main() 
