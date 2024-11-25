#   // INITIAL DESCRIPTION //
"""Run QM Calculations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-09-19'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os
import sys
import subprocess
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from gmx2qmmm.logging import Logger
from gmx2qmmm.generators._helper import filter_xyzq, _flatten
from gmx2qmmm.generators.geometry import read_gmx_structure_header, read_gmx_structure_atoms, read_gmx_box_vectors, write_g96

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class QM():

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

        # XX AJ nma flag needed for some function
        self.nmaflag = 0

    '''
    ------------------------------------
    '''
    # XX Florian and Alina think seperating the functions keeps better flexibility (:
    def generate_qm_input(self):

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
        # XX AJ needed for every qm software
        #   Filter xyzq For Coordinates (No Charges)
        fullcoords = filter_xyzq(self.system.array_xyzq_current, list(range(1,self.system.int_number_atoms)), charges=False)

        atoms_section = ""
        for atom_index, element in enumerate(fullcoords):
            if int(atom_index + 1) in np.array(self.system.list_atoms_qm).astype(int):
                atoms_section += "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                    str(self.system.list_atom_elements[atom_index]),
                    float(element[0]),
                    float(element[1]),
                    float(element[2])
                )
            atom_index += 1

        for element in self.system.list_coordinates_linkatoms:
            atoms_section += "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                str("H"), float(element[0]), float(element[1]), float(element[2])
            )

        atoms_section += "\n"

        # Write the input file
        with open(self.str_inputfile_qm, 'w') as str_inputfile_qm:
            str_inputfile_qm.write(self.header)
            str_inputfile_qm.write(atoms_section)
            str_inputfile_qm.write(self.additional_input)

        if self.dict_input_userparameters['qmcommand'] == 'rung16':



            # XXX Florian!

            # XX AJ I don't understand this stdout stuff, Florian could you check if that's correct?
            original_stdout = sys.stdout

            # XX AJ what is the benefit of using stdout instead of f.write?
            with open('top_info.txt', 'w') as f:
                sys.stdout = f # Change the standard output to the file we created.
                print('QMMMTOPINFO')
                for i in self.system.list_atom_elements:
                    print(str(i))
                sys.stdout = original_stdout # Reset the standard output to its original value
            with open('gro_info.txt', 'w') as f:
                sys.stdout = f # Change the standard output to the file we created.
                print('GROINFO')
                for i in fullcoords:
                    print(str(i))
                sys.stdout = original_stdout # Reset the standard output to its original value

            # Ende Florian


        elif self.dict_input_userparameters['qmcommand'] == 'orca':
            # XX AJ: Nicola can here add the other input files or change the previous part to work with all of the QM programs
            pass

    def run_qm_job(self):

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

        if self.dict_input_userparameters['qmcommand'] == 'g16':
            insert = ""
            if int(self.system.int_step_current) > 0:
                insert = str("." + str(int(self.system.int_step_current)))

            if not os.path.isfile(str(self.str_inputfile_qm) + ".log"):
                # logger(logfile, "Running G16 file.\n")
                self.execute_g16(self.dict_input_userparameters['qmcommand'], str(self.str_inputfile_qm))
                logname = self.str_inputfile_qm[:-3]
                logname += "log"
                os.rename(logname, str(self.dict_input_userparameters['jobname'] + insert + ".gjf.log"))
                os.rename("fort.7", str(self.dict_input_userparameters['jobname'] + insert + ".fort.7"))
                # logger(logfile, "G16 Done.\n")
            else:
                # logger(
                #     logfile,
                #     "NOTE: Using existing G16 files, skipping calculation for this step.\n",
                # )
                pass
            if not os.path.isfile(self.dict_input_userparameters['jobname'] + insert + ".fort.7"):
                if not os.path.isfile("fort.7"):
                    # logger(
                    #     logfile,
                    #     "No fort.7 file was created by the last Gaussian run! Exiting.\n",
                    # )
                    # exit(1)
                    pass
                os.rename("fort.7", str(self.dict_input_userparameters['jobname'] + insert + ".fort.7"))
                # logger(
                #     logfile,
                #     "WARNING: Had to rename fort.7 file but not the log file. MAKE SURE THAT THE FORT.7 FILE FITS TO THE LOG FILE!\n",
                # )
        elif self.dict_input_userparameters['qmcommand'] == 'orca':
            # XX AJ: Nicola adding
            pass

    def execute_g16(self, g16cmd, qmfile):
        #   Call g16 in a seperate function to be able to mock it
        subprocess.call([g16cmd, str(qmfile)])

    def read_qm_energy(self):

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

        if self.dict_input_userparameters['qmcommand'] == 'g16':
            # logger(logfile, "Extracting QM energy.\n")
            self.qmenergy = 0.0
            self.qm_corrdata = []
            if str(self.dict_input_userparameters['program']) == "G16":
                with open(str(self.str_inputfile_qm + ".log")) as ifile:
                    for line in ifile:

                        match = []
                        match2 = []
                        match2 = re.search(
                            r"\sTD[=(\s]", self.dict_input_userparameters['extra'].upper(), flags=re.MULTILINE
                        )
                        if not match2:
                            match2 = re.search(
                                r"^TD[=(\s]", self.dict_input_userparameters['extra'].upper(), flags=re.MULTILINE
                            )
                        if not match2:
                            match2 = re.search(
                                r"\sTD$", self.dict_input_userparameters['extra'].upper(), flags=re.MULTILINE
                            )
                        if not match2:
                            match2 = re.search(
                                r"^TD$", self.dict_input_userparameters['extra'].upper(), flags=re.MULTILINE
                            )
                        if not match2:
                            match = re.search(
                                r"^\s*SCF\s*Done:\s*E\(\S+\)\s*\=\s*([-]*\d+\.\d+)",
                                line,
                                flags=re.MULTILINE,
                            )
                        else:
                            match = re.search(
                                r"^\s*Total\s*Energy,\s*E\(\S+\)\s*\=\s*([-]*\d+\.\d+)",
                                line,
                                flags=re.MULTILINE,
                            )
                        if match:
                            # logger(logfile, "Obtaining charge self-interaction...\n")
                            self.read_pcf_self()
                            # logger(
                            #     logfile, "Done: {:>20.10f} a.u.\n".format(float(pcf_self_pot))
                            # )
                            # G16 energy needs to be corrected for self potential of PCF
                            self.qmenergy = float(match.group(1)) - float(self.pcf_self)
                        match = re.search(r"^\s*ESP\s*charges:", line, flags=re.MULTILINE)
                        if match:
                            for line in ifile:
                                break
                            for line in ifile:
                                match = re.search(
                                    r"^\s*(\d+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                                )
                                if match:
                                    self.qm_corrdata.append(
                                        [
                                            int(match.group(1)),
                                            match.group(2),
                                            float(match.group(3)),
                                        ]
                                    )
                                else:
                                    break
                            break
            # logger(logfile, "QM energy is " + str(float(qmenergy)) + " a.u..\n")
        if self.dict_input_userparameters['qmcommand'] == 'orca':
            # XX AJ: Nicola add
            pass


    def read_qm_forces(self):

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

        self.qmforces = []
        qmonlyforcelist = []
        pcf_grad = []

        if self.dict_input_userparameters['program'] == "G16":
            insert = ""
            if (int(self.system.int_step_current) != 0):
                insert = str("." + str(self.system.int_step_current))
            # logger(logfile,"Reading QM forces using file: "+str(jobname + insert + ".gjf.log")+" and "+ str(jobname + insert + ".fort.7")+"\n")
            qmlogfile = str(self.dict_input_userparameters['jobname'] + insert + ".gjf.log")
            fortfile = str(self.dict_input_userparameters['jobname'] + insert + ".fort.7")

            with open(fortfile) as ifile:
                for line in ifile:
                    match = re.search(
                        r"^\s*(\S*)\s*(\S*)\s*(\S*)", line, flags=re.MULTILINE
                    )
                    if match:
                        qmline = [
                            float(str(match.group(1)).replace("D", "e")) * -1.0,
                            float(str(match.group(2)).replace("D", "e")) * -1.0,
                            float(str(match.group(3)).replace("D", "e")) * -1.0,
                        ]
                        qmonlyforcelist.append(qmline)
            with open(qmlogfile) as i2file:
                for line in i2file:
                    match = re.search(
                        r"^\s*Electrostatic\s*Properties\s*\(Atomic\s*Units\)",
                        line,
                        flags=re.MULTILINE,
                    )
                    if match:
                        for line in i2file:
                            match = re.search(
                                r"^\s*\S+\s*[-]*\d+\.\d+\s*([-]*\d+\.\d+)\s*([-]*\d+\.\d+)\s*([-]*\d+\.\d+)",
                                line,
                                flags=re.MULTILINE,
                            )
                            if match:
                                pcf_grad_line = [
                                    float(match.group(1)),
                                    float(match.group(2)),
                                    float(match.group(3)),
                                ]
                                pcf_grad.append(pcf_grad_line)
                            match = re.search(r"^\s*Leave Link", line, flags=re.MULTILINE)
                            if match:
                                break
                        break
        with open(self.class_topology_qmmm.qmmm_topology) as ifile:
            for line in ifile:
                match = re.search(r"\[\s+moleculetype\s*\]", line)
                if match:
                    for line in ifile:
                        match = re.search(r"\[\s+atoms\s*\]", line)
                        if match:
                            count = 0
                            qmcount = 0
                            m1count = 0
                            for line in ifile:
                                match = re.search(
                                    r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-]*\d+\.*\d*)\s+(\d+\.*\d*)",
                                    line,
                                    flags=re.MULTILINE,
                                )
                                if (
                                    match
                                    and (
                                        int(match.group(1))
                                        not in np.array(self.system.list_atoms_qm).astype(int)
                                    )
                                    and (int(match.group(1)) not in np.array(self.system.list_atoms_m1).astype(int))
                                ):
                                    curr_charge = float(match.group(7))
                                    self.qmforces.append(
                                        [
                                            pcf_grad[count][0] * curr_charge,
                                            pcf_grad[count][1] * curr_charge,
                                            pcf_grad[count][2] * curr_charge,
                                        ]
                                    )
                                    count += 1
                                elif match and int(match.group(1)) in np.array(
                                    self.system.list_atoms_qm
                                ).astype(int):
                                    self.qmforces.append(qmonlyforcelist[qmcount])
                                    qmcount += 1
                                elif match and int(match.group(1)) in np.array(self.system.list_atoms_m1).astype(
                                    int
                                ):
                                    self.qmforces.append(
                                        qmonlyforcelist[m1count + len(self.system.list_atoms_qm)]
                                    )
                                    m1count += 1
                                match = re.search(r"^\s*\n", line, flags=re.MULTILINE)
                                if match:
                                    break
                            break
                    break


    def read_pcf_self(self):

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

        self.pcf_self = 0.0
        with open(self.str_inputfile_qm + ".log") as ifile:
            for line in ifile:
                match = re.search(
                    r"^\s+Self\s+energy\s+of\s+the\s+charges\s+=\s+([-]*\d+\.\d+)\s+a\.u\.",
                    line,
                    flags=re.MULTILINE,
                )
                if match:
                    self.pcf_self = float(match.group(1))
                    break


class QM_gaussian(QM):
    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, str_directory_base):
        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.str_directory_base = str_directory_base

        # XX AJ NMA
        self.nmaflag = 0



    def generate_filename(self):
        insert = ""
        oldinsert = ""
        if int(self.system.int_step_current) > 0:
            insert = str("." + str(int(self.system.int_step_current) ))
            if int(self.system.int_step_current) > 1:
                oldinsert = str("." + str(int(self.system.int_step_current) - 1))
        self.str_inputfile_qm = str(self.dict_input_userparameters['jobname'] + insert + ".gjf")
        self.chkfile = str(self.dict_input_userparameters['jobname'] + insert + ".chk")
        self.oldchkfile = str(self.dict_input_userparameters['jobname'] + oldinsert + ".chk")

    def generate_gaussian_header(self):
        self.header = ""
        self.header += "%NPROCSHARED=" + str(self.dict_input_userparameters['cores']) + "\n"
        self.header += "%MEM=" + str(self.dict_input_userparameters['memory']) + "MB\n"
        self.header += "%CHK=" + self.chkfile + "\n"

        if int(self.system.int_step_current) != 0 or self.nmaflag == 1:
            self.header += "%OLDCHK=" + self.oldchkfile + "\n"

        self.header += "#P " + str(self.dict_input_userparameters['method'])
        self.header += "/" + str(self.dict_input_userparameters['basis'])
        if not self.dict_input_userparameters['extra'].lower() == 'none':
            self.header += " " + str(self.dict_input_userparameters['extra'])

        if int(self.system.int_step_current) != 0 or self.nmaflag == 1:
            self.header += " guess=read"

        self.header += " nosymm gfinput gfprint force charge guess=huckel punch=derivatives iop(3/33=1,6/7=3) prop(field,read) pop=esp\n"
        self.header += "\nQM part of step " + str(self.system.int_step_current) + "\n\n"
        self.header += str(int(self.dict_input_userparameters['charge'])) + " " + str(int(self.dict_input_userparameters['multiplicity'])) + "\n"

    def generate_additional_input(self):
        list_info_pcf = []
        self.additional_input = ""

        with open(self.PCF.pcf_filename) as pcffile:
            for line in pcffile:
                match = re.search(
                    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                )
                if match:
                    # Store the matched values as a list of floats in self.additional_input
                    list_info_pcf.append([
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                        float(match.group(4))
                    ])
        for values in list_info_pcf:
            self.additional_input += "{:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                values[0], values[1], values[2], values[3]
            )

        self.additional_input += '\n'

        for values in list_info_pcf:
            self.additional_input += "{:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                values[0], values[1], values[2]
            )
