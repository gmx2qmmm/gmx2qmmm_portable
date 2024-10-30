#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Run Singlepoint Calculations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-08-12'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os
import math
import sqlite3
import sys
import subprocess
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from Logging.Logger import Logger
from Logging.WriteOutput import Output
from Jobs import QM_job, MM_job
from Generators._helper import filter_xyzq, _flatten
from Generators.GeneratorGeometries import read_gmx_structure_header, read_gmx_structure_atoms, read_gmx_box_vectors, write_g96
from Generators.GeneratorEnergies import GeneratorEnergies, GeneratorForces

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class Singlepoint():

    '''
    This Class Performs A Singlepoint Calculation
    '''

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, str_directory_base) -> None:
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        XX\\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        input_dict: dict -> Dictionary Of The User Parameter Input \\
        class_system: class -> Class Object Of The System \\
        class_topology: class -> Class Object Of The Topology \\
        class_pcf: class -> Class Object Of The Pointchargefield \\
        str_directory_base: str -> Directory Path \\ 
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.str_directory_base = str_directory_base

        #   XX AJ check how to deal with nma flag later
        self.nmaflag = 0

        #   Initialize QM Class (XX AJ differentiate between qm program here?)
        if self.dict_input_userparameters['qmcommand'] == 'g16':
            self.class_qm_job = QM_job.QM_gaussian(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.str_directory_base)
        elif self.dict_input_userparameters['qmcommand'] == 'orca':
            pass

        #   Initialize MM Class
        self.class_mm_job = MM_job.MM(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.str_directory_base)

        #   Initialize QMMM Energy And Forces Generator Classes
        self.class_qmmm_energy = GeneratorEnergies(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.str_directory_base, self.class_qm_job, self.class_mm_job)
        self.class_qmmm_forces = GeneratorForces(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.str_directory_base, self.class_qm_job, self.class_mm_job)


        #   Run (Initial) Singlepoint Calculation
        self.run_calculation()

        #   Generating Output Files
        class_output = Output()
        class_output.oenergy_append(0, self.class_qm_job.qmenergy, self.class_mm_job.mmenergy, self.linkcorrenergy, self.total_energy)
        class_output.oforces_append(0, self.total_force)

    
    def run_calculation(self) -> None:
        
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        XX \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''
        
        # prepare QM input depending on software
        self.class_qm_job.generate_filename()
        self.class_qm_job.generate_gaussian_header()
        self.class_qm_job.generate_additional_input()

        self.class_qm_job.generate_qm_input()

        # run QM calculation
        self.class_qm_job.run_qm_job()

        # read output
        self.class_qm_job.read_qm_energy()  # qm_corrdata correct
        self.class_qm_job.read_qm_forces()

        # prepare MM input
        self.class_mm_job.make_gmx_inp()

        # run MM calculation
        self.class_mm_job.run_gmx()

        # read mm output
        self.class_mm_job.read_mm_energy()
        self.class_mm_job.read_mm_forces()

        #   Calculate Total Energy
        self.linkcorrenergy = self.class_qmmm_energy.get_linkcorrenergy()
        self.total_energy = self.class_qm_job.qmenergy + self.class_mm_job.mmenergy - self.linkcorrenergy

        #   Calculate Total Forces
        self.linkcorrforces = self.class_qmmm_forces.get_linkcorrforces()
        self.total_force = np.array(self.class_qm_job.qmforces) + np.array(self.class_mm_job.mmforces) - np.array(self.linkcorrforces)
        




