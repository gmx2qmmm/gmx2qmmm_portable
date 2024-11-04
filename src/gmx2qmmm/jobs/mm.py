#   // INITIAL DESCRIPTION //
"""Run Gromacs Calculations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-09-19'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from gmx2qmmm.generators.geometry import read_gmx_structure_header, read_gmx_structure_atoms, read_gmx_box_vectors, write_g96

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class MM():

    '''
    This Class Performs A Singlepoint Calculation
    '''

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, str_directory_base) -> None:
        # XX AJ check later which ones of these I actually need here
        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.str_directory_base = str_directory_base

        #   Initialize Gromacs Input
        self.string_structure_gmx_header = read_gmx_structure_header(self.dict_input_userparameters['coordinatefile'])
        self.list_structure_gmx_atoms = read_gmx_structure_atoms(self.dict_input_userparameters['coordinatefile'])
        self.list_box_vectors_initial = read_gmx_box_vectors(self.dict_input_userparameters['coordinatefile'])

        #   Calculate The Maximum Eucledian Distance Between Any Two Atoms And Get New Box Vectors
        array_coordinates_all = self.system.array_xyzq_current[:,:3]
        self.float_distance_max = np.max(np.linalg.norm(array_coordinates_all[np.newaxis, :, :] - array_coordinates_all[:, np.newaxis, :], axis=-1))
        self.list_box_vectors_large = (np.array(self.list_box_vectors_initial) + 10.0 * self.float_distance_max).tolist()



        #   Initialize Gromacs Filenames
        self.mdpname = str(self.dict_input_userparameters['jobname'] + ".mdp")
        self.groname = str(self.dict_input_userparameters['jobname'] + ".boxlarge.g96")
        self.ndxname = str(self.class_topology_qmmm.qmmm_topology + ".ndx")
        self.tprname = str(self.dict_input_userparameters['jobname'] + ".tpr")
        self.trrname = str(self.dict_input_userparameters['jobname'] + ".trr")
        self.xtcname = str(self.dict_input_userparameters['jobname'] + ".xtc")
        self.outname = str(self.dict_input_userparameters['jobname'] + ".out.gro")
        self.gmxlogname = str(self.dict_input_userparameters['jobname'] + ".gmx.log")
        self.edrname = str(self.dict_input_userparameters['jobname'] + ".edr")


    def make_gmx_inp(self):

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

        write_g96(self.groname, self.string_structure_gmx_header, self.list_structure_gmx_atoms, self.system.array_xyzq_current, self.list_box_vectors_large)

        self.prefix =  self.dict_input_userparameters['gmxpath'] + self.dict_input_userparameters['gmxcmd']

        self.write_mdp()

        # XX AJ commented out until testing
        # subprocess.call(
        #     [
        #         prefix,
        #         "grompp",
        #         "-p",
        #         str(self.class_topology_qmmm.qmmm_topology),
        #         "-c",
        #         str(self.groname),
        #         "-n",
        #         str(self.ndxname),
        #         "-f",
        #         str(self.mdpname),
        #         "-o",
        #         str(self.tprname),
        #         "-backup",
        #         "no",
        #     ]
        # )
        # subprocess.call(["rm", "mdout.mdp"])


    def write_mdp(self):

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

        if self.dict_input_userparameters['rcoulomb'] == 0:
            self.dict_input_userparameters['rcoulomb'] = self.float_distance_max
        if self.dict_input_userparameters['rvdw'] == 0:
            self.dict_input_userparameters['rvdw'] = self.dict_input_userparameters['rcoulomb']

        with open(self.mdpname, "w") as ofile:
            ofile.write(
                "title               =  Yo\ncpp                 =  /usr/bin/cpp\nconstraints         =  none\nintegrator          =  md\ndt                  =  0.001 ; ps !\nnsteps              =  1\nnstcomm             =  0\nnstxout             =  1\nnstvout             =  1\nnstfout             =   1\nnstlog              =  1\nnstenergy           =  1\nnstlist             =  1\nns_type             =  grid\nrlist               =  "
            )
            ofile.write(str(float(self.dict_input_userparameters['rcoulomb'])))
            ofile.write(
                "\ncutoff-scheme = group\ncoulombtype    =  cut-off\nrcoulomb            =  "
            )
            ofile.write(str(float(self.dict_input_userparameters['rcoulomb'])))
            ofile.write("\nrvdw                =  ")
            ofile.write(str(float(self.dict_input_userparameters['rvdw'])))
            if self.dict_input_userparameters['useinnerouter']:
                ofile.write(
                "\nTcoupl              =  no\nfreezegrps          =  OUTER\nfreezedim           =  Y Y Y\nenergygrps          =  QM INNER OUTER\nenergygrp-excl = QM QM INNER OUTER OUTER OUTER\nPcoupl              =  no\ngen_vel             =  no\n"
                )
            else:
                ofile.write(
                "\nTcoupl              =  no\nenergygrps          =  QM\nenergygrp-excl = QM QM\nPcoupl              =  no\ngen_vel             =  no\n"
                )

    def run_gmx(self):

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

        pass
        # logger(logfile, "Running Gromacs file.\n")

        # XX AJ commented out until testing
        # subprocess.call(
        # [
        #     self.prefix,
        #     "mdrun",
        #     "-s",
        #     self.tprname,
        #     "-o",
        #     self.trrname,
        #     "-c",
        #     self.outname,
        #     "-x",
        #     self.xtcname,
        #     "-g",
        #     self.gmxlogname,
        #     "-e",
        #     self.edrname,
        #     "-backup",
        #     "no",
        # ]
        # )
        # # os.remove(outname)


    def read_mm_energy(self):

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

        prefix =  self.dict_input_userparameters['gmxpath'] + self.dict_input_userparameters['gmxcmd']

        self.mmenergy = 0.0
        # logger(logfile, "Extracting MM energy.\n")
        # p = subprocess.Popen(
        #     [
        #         prefix,
        #         "energy",
        #         "-f",
        #         edrname,
        #         "-o",
        #         str(edrname + ".xvg"),
        #         "-backup",
        #         "no",
        #     ],
        #     stdout=subprocess.PIPE,
        #     stdin=subprocess.PIPE,
        #     stderr=subprocess.STDOUT,
        # )
        # p.communicate(input=b"11\n\n")

        with open(str(self.edrname + ".xvg")) as ifile:
            for line in ifile:
                match = re.search(
                    r"^    0.000000\s*([-]*\d+.\d+)\n", line, flags=re.MULTILINE
                )
                if match:
                    self.mmenergy = float(match.group(1)) * 0.00038087988
                    break
        # logger(logfile, "MM energy is " + str(float(mmenergy)) + " a.u..\n")




    def read_mm_forces(self):

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

        prefix =  self.dict_input_userparameters['gmxpath'] + self.dict_input_userparameters['gmxcmd']

        self.mmforces = []
        insert = ""
        if int(self.system.int_step_current) != 0:
            insert = str("." + str(self.system.int_step_current))
        # logger(logfile,"Reading MM forces using file: "+str(self.dict_input_userparameters['jobname'] + insert + ".trr/.tpr/.tpr")+"\n")
        trrname = str(self.dict_input_userparameters['jobname'] + insert + ".trr")
        tprname = str(self.dict_input_userparameters['jobname'] + insert + ".tpr")
        xvgname = str(self.dict_input_userparameters['jobname'] + insert + ".xvg")
        # p = subprocess.Popen(
        #     [
        #         prefix,
        #         "traj",
        #         "-fp",
        #         "-f",
        #         trrname,
        #         "-s",
        #         tprname,
        #         "-of",
        #         xvgname,
        #         "-xvg",
        #         "none",
        #         "-backup",
        #         "no",
        #     ],
        #     stdout=subprocess.PIPE,
        #     stdin=subprocess.PIPE,
        #     stderr=subprocess.STDOUT,
        # )
        # p.communicate(input=b"0\n")

        with open(xvgname) as ifile:
            for line in ifile:
                forcelist = re.findall("\S+", line)
                count = 0
                mmforceline = []
                for i in range(1, len(forcelist)):
                    mmforceline.append(float(forcelist[i]) * 2.0155295e-05)
                    count += 1
                    if count > 2:
                        count = 0
                        self.mmforces.append(mmforceline)
                        mmforceline = []
                break  # read only one line
