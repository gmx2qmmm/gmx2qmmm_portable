""" 
Run me with 
python -m unittest -v unit_tests.py
"""
import unittest
import os
import gmx2qmmm.cli
import numpy as np
from unittest.mock import Mock, patch
from gmx2qmmm.generators.geometry import make_xyzq
from gmx2qmmm.jobs.qm import QM
# from gmx2qmmm.jobs.mm import execute_gmx, execute_gmx_communicate

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
    def __init__(self):
        self.correct_output_path = 'ref_output'

    @staticmethod
    def execute_gmx2qmmm_mocking(directory = None):
        if not directory:
# os.chdir(os.path.join('tests', 'opt_input'))
            # os.chdir(os.path.join('tests', 'opt_input'))
            pass
        else:
            os.chdir(directory)

        execute_g16_mock = Mock()
        execute_gmx_mock = Mock()
        execute_gmx_com_mock = Mock()

        with patch('gmx2qmmm.jobs.qm.QM.execute_g16', execute_g16_mock):
            with patch('gmx2qmmm.jobs.mm.MM.run_gmx', execute_gmx_mock):
                with patch('gmx2qmmm.jobs.mm.MM.execute_gmx_communicate', execute_gmx_com_mock):
                    gmx2qmmm.cli.main()
        
        return None
    
    def read_oenergy(self, path=''):
        with open(os.path.join(path, 'oenergy.txt'), 'r') as oenergy:
            lines = oenergy.readlines()
            energies = [i.split('\t') for i in lines[1:]]
            floats = np.array([[float(i.strip('\n')) for i in energies[j]] for j in range(len(energies))])
            values = floats[:, 1:]
            return values
        
    def read_oforces(self, path=''):
        with open(os.path.join(path, 'oforces.txt'), 'r') as oforces:
            lines = oforces.readlines()
            data = [line.split() for line in lines]
            number_of_atoms = data.index(['Step1']) - 2
            data = [line[1:] for line in data]
            data = [row for row in data if row]
            result_array = np.array(data, dtype=float)
            number_of_steps = len(result_array) // number_of_atoms
            columns_per_step = len(result_array[0])
            result_array = result_array.reshape(number_of_steps, number_of_atoms, columns_per_step)
            return result_array

    def compare_oenergy(self):
        calculated_energy = self.read_oenergy()
        correct_energy = self.read_oenergy(path = self.correct_output_path)
        are_equal = np.allclose(calculated_energy, correct_energy)
        self.assertTrue(are_equal)

    def compare_oforces(self):
        calculated_forces = self.read_oforces()
        correct_forces = self.read_oforces(path = self.correct_output_path)
        are_equal = np.allclose(calculated_forces, correct_forces)
        self.assertTrue(are_equal)

    def analyze(self):
        self.compare_oenergy()
        self.compare_oforces()



if __name__ == '__main__':

    testi = TestProductionLevel()
    testi.execute_gmx2qmmm_mocking()
    # unittest.main() 
