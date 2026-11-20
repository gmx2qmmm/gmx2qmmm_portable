#   // INITIAL DESCRIPTION //
"""Run Optimizations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-08-20'

#   // IMPORTS //

#   Imports Of Existing Libraries
import os
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from gmx2qmmm.jobs.singlepoint import Singlepoint
from gmx2qmmm.generators.geometry import propagate_dispvec
from gmx2qmmm.generators._helper import mask_atoms

#   // TODOS & NOTES //
#   TODO:
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class Optimisation():

    '''
    This Class Performs An Optimization
    '''

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, work_dir, base_dir) -> None:
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
        self.work_dir = work_dir
        self.base_dir = base_dir

        self.STEPLIMIT = 0
        self.FTHRESH = 1
        self.STEPSIZE = 2

        #   Perform Initial Singlepoint Calculation
        self.singlepoint = Singlepoint(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.work_dir, self.base_dir)

        #   Setting Up Variables
        self.list_forces_max_all_steps = []
        self.list_forces_all_steps = []
        self.singlepoint.total_force = mask_atoms(self.singlepoint.total_force, self.system.list_atoms_active)
        self.list_forces_all_steps.append(self.singlepoint.total_force)

        self.list_energies_all_steps = []
        self.list_energies_all_steps.append([self.singlepoint.class_qm_job.qmenergy, self.singlepoint.class_mm_job.mmenergy, self.singlepoint.linkcorrenergy, self.singlepoint.total_energy])

        #   Initialize xyzq List, Always Keep The Last Two xyzq In The List -> Therefore, We Start With Two Times The Initial xyzq
        self.list_xyzq_all_steps = [self.system.array_xyzq_current, self.system.array_xyzq_current]

        #   Check If Maximum Force Is Below Threshold
        float_force_max = max(np.max(self.singlepoint.total_force), np.min(self.singlepoint.total_force), key=abs)
        if abs(float_force_max) < float(self.dict_input_userparameters['f_thresh']):
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
        self.PCF.pcf_filename = self.work_dir / (str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".pointcharges"))

        #   Calculate And Apply Displacement
        self.generate_displacement()
        self.system.array_xyzq_current = self.list_xyzq_all_steps[-1]
        self.system.update_linkatom_coordinates()

        #   Update Pointchargefield
        self.PCF.make_pcf()

        #   Run Singlepoint Calculation
        self.singlepoint.run_calculation()

        #    Remove Forces On Inactive Atoms
        self.singlepoint.total_force = mask_atoms(self.singlepoint.total_force, self.system.list_atoms_active)

        #   Evaluate current step
        list_energies_new = [self.singlepoint.class_qm_job.qmenergy, self.singlepoint.class_mm_job.mmenergy, self.singlepoint.linkcorrenergy, self.singlepoint.total_energy]
        self.evaluate_step_steep(list_energies_new)

    def evaluate_step_steep(self, list_energies_new):
        if list_energies_new[-1] > self.list_energies_all_steps[-1][-1]:
            #   Step Gets Rejected, Decrease Stepsize
            self.dict_input_userparameters['stepsize'] = float(self.dict_input_userparameters['stepsize']) * 0.2

            #   Remove Files
            self.remove_previous_files()

            #   Remove Rejected xyzq And Restore Last Accepted Geometry
            self.list_xyzq_all_steps.pop()
            self.system.array_xyzq_current = self.list_xyzq_all_steps[-1]
            self.system.update_linkatom_coordinates()

            #   Rejected Current Step And Use Previous Parameters
            self.system.int_step_current -= 1

            if self.dict_input_userparameters['stepsize'] < 1e-6: #0.000001 a.u.
                self.bool_opt_done = self.STEPSIZE

        else:
            self.dict_input_userparameters['stepsize'] = float(self.dict_input_userparameters['stepsize']) * 1.2

            #   Remove xyzq From Two Steps Ago
            self.list_xyzq_all_steps.pop(0)

            #   Update Forces And Energies
            self.list_forces_all_steps.append(self.singlepoint.total_force)
            self.list_energies_all_steps.append(list_energies_new)

            float_force_max = max(np.max(self.singlepoint.total_force), np.min(self.singlepoint.total_force), key=abs)
            self.list_forces_max_all_steps.append(float_force_max)

            #   Write Output 
            self.singlepoint.class_output.oenergy_append(self.system.int_step_current, self.singlepoint.class_qm_job.qmenergy, self.singlepoint.class_mm_job.mmenergy, self.singlepoint.linkcorrenergy, self.singlepoint.total_energy)
            self.singlepoint.class_output.oforces_append(self.system.int_step_current, self.singlepoint.total_force)

            if abs(float_force_max) < float(self.dict_input_userparameters['f_thresh']):
                self.bool_opt_done = self.FTHRESH

            

    def update_input_filenames(self):
        self.singlepoint.class_mm_job.groname = self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".g96")
        self.singlepoint.class_mm_job.tprname = self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".tpr")
        self.singlepoint.class_mm_job.trrname = self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".trr")
        self.singlepoint.class_mm_job.xtcname = self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".xtc")
        self.singlepoint.class_mm_job.outname = self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".out.gro")
        self.singlepoint.class_mm_job.gmxlogname = self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gmx.log")
        self.singlepoint.class_mm_job.edrname = self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".edr")
        self.singlepoint.class_mm_job.edr_xvg_file = self.work_dir / (str(self.singlepoint.class_mm_job.edrname) + ".xvg")

    def generate_displacement(self):
        if self.dict_input_userparameters['jobtype'] == "SCAN" :
            pass # XX AJ add scan later
            # dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, self.system.int_step_current, True, qmmmInputs.scan_atoms)
        else :
            dispvec = propagate_dispvec(self.dict_input_userparameters['propagator'], self.list_xyzq_all_steps, self.list_forces_all_steps, self.list_forces_max_all_steps[-1], self.dict_input_userparameters['stepsize'], self.system.int_step_current)
        #   Apply Displacement
        list_xyzq_new = self.list_xyzq_all_steps[-1] + np.append(dispvec, np.zeros((len(dispvec),1)), axis=1)
        self.list_xyzq_all_steps.append(list_xyzq_new)

    def remove_previous_files(self):
        os.remove(self.singlepoint.class_mm_job.trrname)
        os.remove(self.singlepoint.class_mm_job.tprname)
        os.remove(self.singlepoint.class_mm_job.gmxlogname)
        os.remove(self.singlepoint.class_mm_job.edrname)
        os.remove(self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".edr.xvg"))
        os.remove(self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gjf.log"))
        os.remove(self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".fort.7"))
        os.remove(self.work_dir / str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gjf"))
        os.remove(self.PCF.pcf_filename)
