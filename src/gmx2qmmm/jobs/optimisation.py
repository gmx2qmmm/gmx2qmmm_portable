#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Run Optimizations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-08-20'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os
import math
import sqlite3
import sys
import copy
import subprocess
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries
import Jobs.Singlepoint as SP

#   Imports From Custom Libraries
from Logging.Logger import Logger
from Generators.GeneratorGeometries import propagate_dispvec
from Generators._helper import filter_xyzq, _flatten, mask_atoms

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class Optimisation():

    '''
    This Class Performs An Optimization
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

        #define done XX AJ evaluate this later (why are we not using just the strings?)
        self.STEPLIMIT = 0
        self.FTHRESH = 1
        self.STEPSIZE = 2

        #   Perform Initial Singlepoint Calculation
        self.singlepoint = SP.Singlepoint(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.str_directory_base)

        #   Setting Up Variables
        self.list_forces_max_all_steps = []
        self.list_forces_all_steps = []
        self.list_forces_all_steps.append(self.singlepoint.total_force)

        self.list_energies_all_steps = []
        self.list_energies_all_steps.append([self.singlepoint.class_qm_job.qmenergy, self.singlepoint.class_mm_job.mmenergy, self.singlepoint.linkcorrenergy, self.singlepoint.total_energy])

        #   Initialize xyzq List, Always Keep The Last Two xyzq In The List -> Therefore, We Start With Two Times The Initial xyzq
        self.list_xyzq_all_steps = [self.system.array_xyzq_current, self.system.array_xyzq_current]

        # XX AJ at this point I don't think we still need this function, but I'll keep these comments for now in case I'm wrong
        # self.force_clean = self.make_clean_force()

        #   Check If Maximum Force Is Below Threshold
        float_force_max = max(np.max(self.singlepoint.total_force), np.min(self.singlepoint.total_force), key=abs)
        if abs(float_force_max) < self.dict_input_userparameters['f_thresh']:
            self.bool_opt_done = self.FTHRESH
        else:
            self.bool_opt_done = False
        self.list_forces_max_all_steps.append(float_force_max)

        #   Starting Optimization Cycles
        while not self.bool_opt_done and self.system.int_step_current <= self.dict_input_userparameters['maxcycle']:
            self.run_optimization_cycle()


    def run_optimization_cycle(self):

        self.system.int_step_current += 1

        #   Prepare New Input
        self.update_input_filenames()

        #   Update Pointchargefield Filename
        self.PCF.pcf_filename = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".pointcharges")

        #   Calculate And Apply Displacement
        self.generate_displacement()
 
        #   Update Pointchargefield
        self.PCF.make_pcf()

        #   Run Singlepoint Calculation
        self.singlepoint.run_calculation()

        #    Remove Forces On Inactive Atoms
        self.singlepoint.total_force = mask_atoms(self.singlepoint.total_force, self.system.list_atoms_active)

        #   Update Forces And Energies
        self.list_forces_all_steps.append(self.singlepoint.total_force)
        self.list_energies_all_steps.append([self.singlepoint.class_qm_job.qmenergy, self.singlepoint.class_mm_job.mmenergy, self.singlepoint.linkcorrenergy, self.singlepoint.total_energy])

        #   The following section needs to be adapted for different optimizers, for bfgs we accept increasing energies, for steep and conjugate we dont but we then also dont need to check for the max force
        #   XX AJ the following code works for steepest descent. I will change it later for different optimizers
        #   Check If The Total Energy Improved

        #   Evaluate current step
        if self.dict_input_userparameters['propagator'] == 'steep':
            self.evaluate_step_steep()
        elif self.dict_input_userparameters['propagator'] == 'conjgrad':
            self.evaluate_step_conjgrad()
        elif self.dict_input_userparameters['propagator'] == 'bfgs':
            self.evaluate_step_bfgs()

        # XX AJ take care of the output later
        if self.dict_input_userparameters['jobname'] == "SCAN" :
            # write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step, energy_file="oenergy_%s.txt"%jobname, forces_file="oforces_%s.txt"%jobname)
            pass
        else:
            # write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step)
            pass
        # gro = qmmmInputs.gro            #SIMON
        # logger(logfile, "Due to the decrease of the energy, the structure "+str(gro)+" will be used from now on.\n")

        
        pass

    def evaluate_step_steep(self):
        if self.list_energies_all_steps[-1][-1] > self.list_energies_all_steps[-2][-1]:
            #   Step Gets Rejected, Decrease Stepsize
            self.dict_input_userparameters['stepsize'] *= 0.2

            #   Remove Files
            self.remove_previous_files()

            #   Rejected Current Step And Use Previous Parameters
            self.system.int_step_current -= 1

            self.list_xyzq_all_steps.pop()
            self.list_forces_all_steps.pop()
            self.list_forces_max_all_steps.pop()

        else:
            self.dict_input_userparameters['stepsize'] *= 1.2

            #   Remove xyzq From Two Steps Ago
            self.list_xyzq_all_steps.pop(0)

            float_force_max = max(np.max(self.singlepoint.total_force), np.min(self.singlepoint.total_force), key=abs)
            self.list_forces_max_all_steps.append(float_force_max)

            if abs(float_force_max) < float(self.dict_input_userparameters['f_thresh']):
                # logger(logfile,"Max force (%f) below threshold (%f) Finishing.\n"%(maxforce,f_thresh))
                self.bool_opt_done = self.FTHRESH

            elif float(self.dict_input_userparameters['stepsize']) < 1e-6: #0.000001 a.u.
                self.bool_opt_done = self.STEPSIZE
                # logger(
                #     logfile,
                #         ("Step became lower than 0.000001 a.u., optimization is considered done for now. " +
                #          "This is the best we can do unless reaching unacceptable numerical noise levels.\n"),
                # )
    
    def evaluate_step_bfgs(self):
        #   Check If Optimization Is Finished
        if abs(self.float_force_max) < float(self.dict_input_userparameters['f_thresh']):
            self.bool_opt_done = self.FTHRESH
        elif float(self.dict_input_userparameters['stepsize']) < 1e-6: #0.000001 a.u.
            self.bool_opt_done = self.STEPSIZE

 
    def update_input_filenames(self):
        self.singlepoint.groname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".g96")
        self.singlepoint.tprname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".tpr")
        self.singlepoint.trrname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".trr")
        self.singlepoint.xtcname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".xtc")
        self.singlepoint.outname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".out.gro")
        self.singlepoint.gmxlogname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gmx.log")
        self.singlepoint.edrname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".edr")

    def generate_displacement(self):
        if self.dict_input_userparameters['jobtype'] == "SCAN" :
            pass # XX AJ add scan later
            # dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, self.system.int_step_current, True, qmmmInputs.scan_atoms)
        else :
            dispvec = propagate_dispvec(self.dict_input_userparameters['propagator'], self.list_xyzq_all_steps, self.list_forces_all_steps, self.list_forces_max_all_steps[-1], self.dict_input_userparameters['stepsize'], self.system.int_step_current)
        #    write_dispvec(dispvec, curr_step, count_trash) Simon implemented this to check if the dispvec is correct!

        #   Apply Displacement
        list_xyzq_new = self.list_xyzq_all_steps[-1] + np.append(dispvec, np.zeros((len(dispvec),1)), axis=1)
        self.list_xyzq_all_steps.append(list_xyzq_new)

    def remove_previous_files(self):
        os.remove(self.singlepoint.trrname)
        os.remove(self.singlepoint.tprname)
        os.remove(self.singlepoint.gmxlogname)
        os.remove(self.singlepoint.edrname)
        os.remove(str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".edr.xvg"))
        os.remove(str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gjf.log"))
        os.remove(str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".fort.7"))
        os.remove(str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gjf"))
        os.remove(self.PCF.pcf_filename)
