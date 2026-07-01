#   // INITIAL DESCRIPTION //
"""Run Singlepoint Calculations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-08-12'

from collections.abc import Mapping
from typing import Any

import numpy as np

from gmx2qmmm.generators.system import SystemInfo
from gmx2qmmm.generators.topology import GenerateTopology
from gmx2qmmm.generators.pcf.base import PCFGenerator
from gmx2qmmm.generators._helper import filter_xyzq, _flatten
from gmx2qmmm.generators.geometry import read_gmx_structure_header, read_gmx_structure_atoms, read_gmx_box_vectors, write_g96
from gmx2qmmm.generators.energies import GeneratorEnergies, GeneratorForces
from gmx2qmmm.logging_utils import Output
from gmx2qmmm.jobs import mm, qm
from gmx2qmmm.types import StrPath

#   // TODOS & NOTES //
#   TODO:
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class Singlepoint():

    '''
    This Class Performs A Singlepoint Calculation
    '''

    def __init__(
        self,
        parameters: Mapping[str, Any],
        system: SystemInfo,
        topology: GenerateTopology,
        *,
        pcf_generator: PCFGenerator,
        work_dir: StrPath,
        base_dir: StrPath
    ) -> None:
        """Initialise Singlepoint Job

        Args:
            parameters: A dictionary containing all input parameters for the job
            system: An instance of :class:`~gmx2qmmm.generators.system.SystemInfo`
                containing information about the system
            topology: An instance of :class:`~gmx2qmmm.generators.topology.GenerateTopology`
            pcf_generator: A point charge field generator
            work_dir: The working directory for the job
            base_dir: The base directory for the job
        """

        self.system = system
        self.topology = topology
        self.parameters = parameters
        self.pcf_generator = pcf_generator
        self.work_dir = work_dir
        self.base_dir = base_dir

        #   TODO: check how to deal with nma flag later
        self.nmaflag = 0

        #   Initialize QM Class (differentiate between qm program here?)
        if self.parameters['qmcommand'] == 'g16':
            self.class_qm_job = qm.QM_gaussian(self.parameters, self.system, self.topology, self.pcf_generator, self.work_dir)
        elif self.parameters['qmcommand'] == 'orca':
            pass

        #   Initialize MM Class
        self.class_mm_job = mm.MM(self.parameters, self.system, self.topology, self.pcf_generator, self.work_dir)

        #   Initialize QMMM Energy And Forces Generator Classes
        self.class_qmmm_energy = GeneratorEnergies(self.parameters, self.system, self.topology, self.pcf_generator, self.work_dir, self.base_dir, self.class_qm_job, self.class_mm_job)
        self.class_qmmm_forces = GeneratorForces(self.parameters, self.system, self.topology, self.pcf_generator, self.work_dir, self.base_dir, self.class_qm_job, self.class_mm_job)


        #   Run (Initial) Singlepoint Calculation
        self.run_calculation()

        #   Generating Output Files
        class_output = Output(self.work_dir)
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
