#   // INITIAL DESCRIPTION //
# NOTE (AJ): I don't know what to put here because the module only has this one class. What is the difference between the module description and the class description here?
"""Short Module Description; Reference To Readme"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-03-25'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os
import math
import json
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from gmx2qmmm.logging import Logger
from gmx2qmmm.generators import geometry, pcf_from_top
from gmx2qmmm.generators._helper import _flatten, get_coordinates_linkatoms_angstrom

#   // TODOS & NOTES //
#   TODO:
#   - Add the option to read atom input from a list (if that's what we want, wait for Florians opinion)
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class SystemInfo():

    '''
    This Class Reads And Stores Information About The System
    '''

    def __init__(self, dict_input_userparameters) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Calls Functions To Read Information About The System\\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        input_dict: dict -> Dictionary Of The User Parameter Input \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        self.dict_input_userparameters = dict_input_userparameters
        self.int_step_current = 0

        #   Make A List Of All Topology Files
        self.list_topology = self.get_list_topologies(self.dict_input_userparameters['topologyfile'])

        # Read The Different Types Of Molecules In The System
        self.list_molecules = self.read_molecules()

        # Read The Charges For All Atoms
        self.list_charges = []
        for molecule in self.list_molecules:
            self.list_charges.extend(self.readcharges(molecule))

        #   Read All Atom Lists
        #   XX AJ currently I'm assuming we're always having separate files for atom indices, only remove this comment when we finally decided so or add the possibility for lists
        self.list_atoms_qm = self.read_atoms_list(self.dict_input_userparameters['qmatomslist'])
        if self.dict_input_userparameters['jobtype'] != 'SINGLEPOINT':
            self.list_atoms_active = self.read_atoms_list(self.dict_input_userparameters['activeatomslist'])
        else:
            self.list_atoms_active = []
        if dict_input_userparameters['useinnerouter']:
            self.list_atoms_inner = self.read_atoms_list(self.dict_input_userparameters['inneratomslist'])
            self.list_atoms_outer = self.read_atoms_list(self.dict_input_userparameters['outeratomslist'])
        #   If we're not using inner/outer we're keeping the lists empty (currently necessary for some functions but that might change later)
        else:
            self.list_atoms_inner = []
            self.list_atoms_outer = []

        #   Read Connectivity
        self.list_connectivity_topology = self.read_connectity_topology(self.dict_input_userparameters['topologyfile'])

        #   Read Initial Geometry
        #   XX AJ I would prefer only one function here independent of the file type and only make that distinction within the function. I'll get back to that later when I'm writing GeneratorGeometries.py
        if self.dict_input_userparameters['coordinatefile'][-4:] == ".gro":
            #logger(logfile, "Reading geometry (.gro)...\n")
            # XX AJ if we rewrite gro files to g96 files, the following function can be deleted and we can read the g96 file with geometry.readg96 afterwards
            self.list_geometry_initial = geometry.readgeo(self.dict_input_userparameters['coordinatefile'])
            # Writing high-precision coordinate file
            # logger(logfile, "Writing high-precision coordinate file...")
            self.write_file_gromacs_highprec(self.dict_input_userparameters['coordinatefile']) # XX temp removed logfile until logfile decision of Florian AJ
            self.dict_input_userparameters['coordinatefile'] = self.dict_input_userparameters['jobname'] + ".g96"
            # logger(logfile, "Done.\n")
        elif self.dict_input_userparameters['coordinatefile'][-4:] == ".g96":
            #logger(logfile, "Reading geometry (.g96)...\n")
            self.list_geometry_initial = geometry.readg96(self.dict_input_userparameters['coordinatefile'])

        self.int_number_atoms = int(len(self.list_geometry_initial)/3)

        #   Create xyzq (Coordinates And Charges) For The Whole System
        #   XX AJ also prefer only one function here, I'll make one function of it
        #   XX AJ technically, I would prefer this xyzq function not in this class, but it's used in the get_linkatoms_ang, so I'll keep it here
        if dict_input_userparameters['useinnerouter']:
            self.array_xyzq_initial = geometry.make_xyzq_io(self.list_geometry_initial, self.list_charges, self.list_atoms_outer)
        else:
            self.array_xyzq_initial = geometry.make_xyzq(self.list_geometry_initial, self.list_charges)

        #   Current Geometry Gets Updated During Optimizations
        self.array_xyzq_current = self.array_xyzq_initial

        #   Read linkatoms and next order atoms in mm region
        self.list_atoms_m1 = self.get_list_atoms_m1()
        self.list_atoms_m2 = self.get_list_atoms_m2()

        #   Read coordinates of linkatoms in angstrom
        self.list_coordinates_linkatoms = get_coordinates_linkatoms_angstrom(self.array_xyzq_initial, self.list_atoms_qm, self.list_atoms_m1, self.list_connectivity_topology, [])
        self.linkcorrlist, self.list_atoms_q1, self.list_atoms_q2, self.list_atoms_q3, self.list_atoms_m3 = self.get_list_atoms_link()



    def read_atoms_list(self, file_input_atoms) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        ---------------
        Reads Atom Indices
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        file_input_atoms: str -> Index File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_atoms: list -> List Of Atom Indices \\
        ------------------------------ \\
        '''

        list_atoms = []
        with open(file_input_atoms, 'r') as atom_file:
            for line in atom_file:

                # Skip Index Group Name
                if "[" in line or "]" in line:
                    continue

                atom_index_list = re.findall("\d+", line)
                if atom_index_list:
                    list_atoms.extend(map(int, atom_index_list))
        list_atoms = sorted(np.array(list_atoms).astype(int))
        return list_atoms

    def read_molecules(self) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads List Of Molecules From The Topology File \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_molecule: list -> List Of Molecules And Their Amount In The System \\
        ------------------------------ \\
        '''

        list_molecule = []
        with open(self.dict_input_userparameters['topologyfile'], 'r') as file_input:
            bool_match = False
            for line in file_input:
                match = re.search(r"^\[ molecules \]", line, flags=re.MULTILINE)
                if match:
                    bool_match = True
                    break
            if not bool_match:
                #   XX AJ turn into Logging
                print('No "molecules" entry in ' + str(self.top) + " found. Exiting.")
                exit(1)
            for line in file_input:
                # Skip All Comment Lines
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                # Find All Molecule Types And Their Amount
                else:
                    # Extract Lines With A String (Molecule Type, Group 1) And A Number (Molecule Amount, Group 2) Seperated By Whitespace
                    match = re.search(r"^(\S+)\s+(\d+)", line, flags=re.MULTILINE)
                    if match:
                        list_molecule.append([match.group(1), match.group(2)])
                    else:
                        #   XX AJ turn into Logging
                        print("Found an incomprehensible line in molecules list. Exiting.")
                        print("Last line was:")
                        print(line)
                        exit(1)
        return list_molecule

    def readcharges(self, list_molecule_entry) -> list:
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads Charges For All Atoms In A Molecule Type \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        list_molecule_entry: list -> List With Molecule Name And Amount Of This Molecule \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_charges: list -> List Of Charges For Every Atom In This Molecule Types \\
        ------------------------------ \\
        '''

        list_charges = []
        current_topology = self.dict_input_userparameters['topologyfile']
        molecule_name = list_molecule_entry[0]
        molecule_count = int(list_molecule_entry[1])
        found = self.check_occurence_topology_molecule(molecule_name, self.dict_input_userparameters['topologyfile'])

        if not found:
            for element in self.list_topology:
                found = self.check_occurence_topology_molecule(molecule_name, element)
                if found:
                    current_topology = element
                    break
        if not found:
            print("No charges found for " + str(molecule_name) + ". Exiting.")
            exit(1)
        with open(current_topology) as ifile:
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s*moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        matchstring = r"^\s*" + re.escape(molecule_name)
                        match = re.search(matchstring, line, flags=re.MULTILINE)
                        if match:
                            found = True

                            for line in ifile:
                                match = re.search(
                                    r"^\[\s*atoms\s*\]", line, flags=re.MULTILINE
                                )
                                if match:
                                    break
                            break
                        else:
                            found = False
                            break
                    if found:
                        break
            for line in ifile:
                match = re.search(r"^\[", line, flags=re.MULTILINE)
                if match:
                    break
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(
                    r"^\s*\d+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+([-]*\d+[\.]*[\d+]*)*",
                    line,
                    flags=re.MULTILINE,
                )
                if match:
                    list_charges.append(float(match.group(1)))

        list_charges *= molecule_count

        return list_charges

    def check_occurence_topology_molecule(self, molecule_name, topology_file) -> bool:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Checks If The Molecule (molecule_name) Is Listed In The Topology File \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        bool_occurence_molecule: bool \\
        ------------------------------ \\
        '''

        with open(topology_file) as ifile:
            bool_occurence_molecule = False
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s*moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        if line.strip() == '' or line.startswith(';'):
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                bool_occurence_molecule = True
                            break
                if bool_occurence_molecule:
                    break

        return bool_occurence_molecule

    def read_connectity_topology(self, topology_file) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of All Molecules In The System \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        connectivity_list: list -> List Of The Connectivity Of All Molecules \\
        ------------------------------ \\
        '''

        int_count = 0
        connectivity_list = []
        for molecule in self.list_molecules:
            current_topology = self.get_current_topology(molecule[0], topology_file)
            for _ in range(0, int(molecule[1])):
                mollength = self.get_mollength_direct(molecule[0], current_topology)
                connset = self.get_list_connectivity(int_count, molecule[0], current_topology)
                connectivity_list += connset
                int_count += int(mollength)

        return connectivity_list

    def get_current_topology(self, molecule_name, topology_file) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of All Molecules In The System \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        connectivity_list: list -> List Of The Connectivity Of All Molecules \\
        ------------------------------ \\
        '''

        str_current_topology = topology_file
        bool_occurence_molecule = pcf_from_top.checkformol(molecule_name, topology_file)

        if not bool_occurence_molecule:
            for element in self.list_topology:
                bool_occurence_molecule = pcf_from_top.checkformol(molecule_name, element)
                if bool_occurence_molecule:
                    str_current_topology = element
                    break
        if not bool_occurence_molecule:
            # XX AJ replace with logger
            print("No charges found for " + str(molecule_name) + ". Exiting.")
            exit(1)

        return str_current_topology

    def get_mollength_direct(self, molecule_name, topology_file) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Amount Of Atoms Of The Molecule (molecule_name) \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        connectivity_list: list -> List Of The Connectivity Of All Molecules \\
        ------------------------------ \\
        '''

        mollength = 0
        with open(topology_file) as ifile:
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                for line in ifile:
                                    match = re.search(r"^;", line, flags=re.MULTILINE)
                                    if match:
                                        continue
                                    else:
                                        match = re.search(
                                            r"^\[\s*atoms\s*\]", line, flags=re.MULTILINE
                                        )
                                        if match:
                                            for line in ifile:
                                                match = re.search(
                                                    r"^;", line, flags=re.MULTILINE
                                                )
                                                if match:
                                                    continue
                                                else:
                                                    match = re.search(
                                                        r"^\s*(\d+)\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+[-]*\d+[\.]*[\d+]*",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if match:
                                                        mollength = int(match.group(1))
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                            break
                                break
                    break
        return mollength

    def get_list_connectivity(self, offset, molecule_name, topology_file) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of One Molecule (molecule_name) \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        offset: XX ?
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        connectivity_list: list -> List Of The Connectivity Of One Molecule \\
        ------------------------------ \\
        '''

        connectivity_list = []
        with open(topology_file) as ifile:
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[ moleculetype \]", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                for line in ifile:
                                    match = re.search(r"^;", line, flags=re.MULTILINE)
                                    if match:
                                        continue
                                    else:
                                        match = re.search(
                                            r"^\[ bonds \]", line, flags=re.MULTILINE
                                        )
                                        if match:
                                            for line in ifile:
                                                match = re.search(
                                                    r"^;", line, flags=re.MULTILINE
                                                )
                                                if match:
                                                    continue
                                                else:
                                                    match = re.search(
                                                        r"^\s*(\d+)\s+(\d+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if match:
                                                        a = int(match.group(1)) + int(
                                                            offset
                                                        )
                                                        b = int(match.group(2)) + int(
                                                            offset
                                                        )
                                                        if a > b:
                                                            c = b
                                                            b = a
                                                            a = c
                                                        found = False
                                                        for element in connectivity_list:
                                                            if int(element[0]) == a:
                                                                if int(b) not in np.array(
                                                                    element
                                                                ).astype(int):
                                                                    element.append(int(b))
                                                                    found = True
                                                        if not found:
                                                            connectivity_list.append(
                                                                [int(a), int(b)]
                                                            )
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                            break
                                break
                    break
        return connectivity_list

    def write_file_gromacs_highprec(self, gro_file) -> str:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Writes A g96 File From The gro File \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        gro_file: str -> Name Of The gro File \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        filename: str -> Name Of The New g96 File \\
        ------------------------------ \\
        '''

        #  XX not tested yet, bc I didn't use a gro file AJ
        filename=str(self.dict_input_userparameters['jobname'] + ".g96")
        with open(filename,"w") as ofile:
                with open(gro_file) as ifile:
                        count=0
                        for line in ifile:
                                count+=1
                                if count==2:
                                        break
                        ofile.write("TITLE\nProtein\nEND\nPOSITION\n")
                        counter=0
                        finalline=""
                        for line in ifile:
                                match=re.search(r'^([\s\d]{5})(.{5})(.{5})([\s\d]{5})\s*([-]*\d+\.\d*)\s*([-]*\d+\.\d*)\s*([-]*\d+\.\d*)', line,flags=re.MULTILINE)
                                if not match:
                                        # logger(logfile,str("Successfully wrote " +  str(int(counter)) + " atoms to internal precision format file.\n"))
                                        finalline=line
                                        break
                                else:
                                        counter+=1
                                        ofile.write(str(match.group(1))+" "+ str(match.group(2))+" "+str(match.group(3))+" {:>6d} ".format(int(counter))+ "{:>15.9f} {:>15.9f} {:>15.9f}\n".format(float(match.group(5)),float(match.group(6)),float(match.group(7))))
                        ofile.write("END\nBOX\n")
                        match=re.search(r'^\s*(\d+\.*\d*)\s*(\d+\.*\d*)\s*(\d+\.*\d*)', finalline,flags=re.MULTILINE)
                        if not match:
                                # logger(logfile,str("Unexpected line instead of box vectors. Exiting. Last line:\n"))
                                # logger(logfile,line)
                                exit(1)
                        else:
                                ofile.write(str(" {:>15.9f} {:>15.9f} {:>15.9f}\n".format(float(match.group(1)),float(match.group(2)),float(match.group(3)))))
                        ofile.write("END")

        return filename

    def get_list_atoms_m1(self) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All M1 Atoms \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        m1list: list -> List Of M1 Atom Indices \\
        ------------------------------ \\
        '''

        m1list = []
        for element in self.list_atoms_qm:
            bondlist = self.get_list_atoms_m1qm(element)
            for entry in bondlist:
                if (int(entry) not in np.array(self.list_atoms_qm).astype(int)) and (
                    int(entry) not in np.array(m1list).astype(int)
                ):
                    m1list.append(entry)
        return m1list

    def get_list_atoms_m1qm(self, target) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All M1 Atoms For A Specific QM Atom \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        target: int -> Index Of Specific QM Atom \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        partnerlist: list -> List Of M1 Atoms For This Specific QM Atom \\
        ------------------------------ \\
        '''

        partnerlist = []
        for entry in self.list_connectivity_topology:
            found = False

            for i in range(0, len(entry)):
                if int(entry[i]) == int(target):
                    found = True
                    break
            if found:
                if int(entry[0]) == int(target):
                    for i in range(1, len(entry)):
                        if int(entry[i]) not in np.array(partnerlist).astype(int):
                            partnerlist.append(int(entry[i]))
                else:
                    if int(entry[0]) not in np.array(partnerlist).astype(int):
                        partnerlist.append(int(entry[0]))

        return partnerlist

    def get_list_atoms_m2(self) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All M2 Atoms \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        m2list: list -> List Of M2 Atom Indices \\
        ------------------------------ \\
        '''

        list_atoms_m2 = []
        for element in self.list_atoms_m1:
            m2line = []
            bondlist = self.get_list_atoms_m1qm(element)
            for entry in bondlist:
                if int(entry) not in np.array(self.list_atoms_qm).astype(int):
                    m2line.append(entry)
            list_atoms_m2.append(m2line)
        return list_atoms_m2

    def get_list_atoms_link(self) -> tuple:
        # XX AJ evaluate return type later with Florian
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All Q1, Q2, Q3 And M3 Atom Indices \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        final_linkcorrlist: list -> List Of All QM-MM-Atom Combinations That Are Excluded For Coulomb Forces \\
        q1list: list -> List Of Q1 Atom Indices \\
        q2list: list -> List Of Q2 Atom Indices \\
        q3list: list -> List Of Q3 Atom Indices \\
        m3list: list -> List Of M3 Atom Indices \\
        ------------------------------ \\
        '''
        linkcorrlist = []
        m3list = []
        m4list = []
        q1list = []
        q2list = []
        q3list = []
        # get q1 and m2
        for element in self.list_atoms_m1:
            q1line = []
            for entry in self.list_connectivity_topology:
                if int(element) in np.array(entry).astype(int):
                    if int(element) != int(entry[0]):
                        if int(entry[0]) in np.array(self.list_atoms_qm).astype(int):
                            q1line.append(int(entry[0]))
                    else:
                        for i in range(1, len(entry)):
                            if int(entry[i]) in np.array(self.list_atoms_qm).astype(int):
                                q1line.append(int(entry[i]))
            q1list.append(q1line)
        # get q2
        q1list = list(_flatten(q1list))
        for element in q1list:
            q2line = []
            for conn in self.list_connectivity_topology:
                if int(element) in np.array(conn).astype(int):
                    if (
                        int(element) != int(conn[0])
                        and (int(conn[0]) in np.array(self.list_atoms_qm).astype(int))
                        and (int(conn[0]) not in np.array([int(x) for x in _flatten(q1list)]))
                    ):
                        q2line.append(int(conn[0]))
                    elif int(element) == int(conn[0]):
                        for i in range(1, len(conn)):
                            if (int(conn[i]) in np.array(self.list_atoms_qm).astype(int)) and (
                                int(conn[i]) not in np.array([int(x) for x in _flatten(q1list)])
                            ):
                                q2line.append(int(conn[i]))
            q2list.append(q2line)
        # get q3
        for element in q2list:
            q3lineline = []
            for entry in element:
                q3line = []
                for conn in self.list_connectivity_topology:
                    if int(entry) in np.array(conn).astype(int):
                        if (
                            int(entry) != int(conn[0])
                            and (int(conn[0]) in np.array(self.list_atoms_qm).astype(int))
                            and int(conn[0]) not in np.array(q1list).astype(int)
                            and int(conn[0]) not in np.array([int(x) for x in _flatten(q1list)])
                        ):
                            q3line.append(int(conn[0]))
                        elif int(entry) == int(conn[0]):
                            for i in range(1, len(conn)):
                                if (
                                    int(conn[i]) in np.array(self.list_atoms_qm).astype(int)
                                    and int(conn[i]) not in np.array(q1list).astype(int)
                                    and int(conn[i]) not in np.array([int(x) for x in _flatten(q1list)])
                                ):
                                    q3line.append(int(conn[i]))
                q3lineline.append(q3line)
            q3list.append(q3lineline)
        # get m3
        for element in self.list_atoms_m2:
            m3lineline = []
            for entry in element:
                m3line = []
                for conn in self.list_connectivity_topology:
                    if int(entry) in np.array(conn).astype(int):
                        if (
                            int(entry) != int(conn[0])
                            and (int(conn[0]) not in np.array(self.list_atoms_qm).astype(int))
                            and int(conn[0]) not in np.array(self.list_atoms_m1).astype(int)
                        ):
                            m3line.append(int(conn[0]))
                        elif int(entry) == int(conn[0]):
                            for i in range(1, len(conn)):
                                if int(conn[i]) not in np.array(self.list_atoms_qm).astype(
                                    int
                                ) and int(conn[i]) not in np.array(self.list_atoms_m1).astype(int):
                                    m3line.append(int(conn[i]))
                m3lineline.append(m3line)
            m3list.append(m3lineline)
        # get m4
        for element in m3list:
            m4linelineline = []
            for entry in element:
                m4lineline = []
                for stuff in entry:
                    m4line = []
                    for conn in self.list_connectivity_topology:
                        if int(stuff) in np.array(conn).astype(int):
                            if (
                                int(stuff) != int(conn[0])
                                and (int(conn[0]) not in np.array(self.list_atoms_qm).astype(int))
                                and int(conn[0]) not in np.array(self.list_atoms_m1).astype(int)
                            ):
                                found = False
                                for morestuff in self.list_atoms_m2:
                                    if int(conn[0]) in np.array(morestuff).astype(int):
                                        found = True
                                        break
                                if not found:
                                    m4line.append(int(conn[0]))
                            elif int(stuff) == int(conn[0]):
                                for i in range(1, len(conn)):
                                    if int(conn[i]) not in np.array(self.list_atoms_qm).astype(
                                        int
                                    ) and int(conn[i]) not in np.array(self.list_atoms_m1).astype(int):
                                        found = False
                                        for morestuff in self.list_atoms_m2:
                                            if int(conn[i]) in np.array(morestuff).astype(
                                                int
                                            ):
                                                found = True
                                                break
                                        if not found:
                                            m4line.append(int(conn[i]))
                    m4lineline.append(m4line)
                m4linelineline.append(m4lineline)
            m4list.append(m4linelineline)
        # set up link atom charge corr pairs: q1-m2, q1-m3, q2-m2, l-m2, l-m3, l-m4 - l are represented as their m1 counterparts!
        # for all these combinations we later have to calculate the coulomb potential and subtract that because they all are 1-4 excluded
        count = 0
        for element in self.list_atoms_m1:
            linkpairline = []
            for entry in self.list_atoms_m2[count]:
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            for entry in _flatten(m3list[count]):
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            for entry in _flatten(m4list[count]):
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            linkcorrlist.append(linkpairline)
            count += 1
        count = 0
        for element in q1list:
            linkpairline = []
            for stuff in self.list_atoms_m2[count]:
                if int(element) < int(stuff):
                    linkpairline.append([element, stuff])
                else:
                    linkpairline.append([stuff, element])
            for stuff in list(_flatten(m3list[count])):
                if int(element) < int(stuff):
                    linkpairline.append([element, stuff])
                else:
                    linkpairline.append([stuff, element])
            linkcorrlist.append(linkpairline)
            count += 1
        count = 0
        for element in q2list:
            for entry in element:
                linkpairline = []
                for stuff in self.list_atoms_m2[count]:
                    if int(entry) < int(stuff):
                        linkpairline.append([entry, stuff])
                    else:
                        linkpairline.append([stuff, entry])
                linkcorrlist.append(linkpairline)
            count += 1
        reshaped_linkcorrlist = np.array(list(_flatten(linkcorrlist))).reshape(-1, 2)
        final_linkcorrlist = []
        for element in reshaped_linkcorrlist:
            found = False
            for i in range(0, len(final_linkcorrlist)):
                if (
                    final_linkcorrlist[i][0] == element[0]
                    and final_linkcorrlist[i][1] == element[1]
                ):
                    found = True
                    break
            if found:
                continue
            final_linkcorrlist.append(element)
        return (
            sorted(final_linkcorrlist, key=lambda l: l[0]),
            q1list,
            q2list,
            q3list,
            m3list,
        )

    def get_list_topologies(self, topology_input) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        ---------------  \\
        Recursively Looks For Other Topology Files Referred To In The Original File \\
        And Adds Them To The List Of Topology Files \\
        ------------------------------ \\
        INPUT: \\
        ---------------  \\
        topology_input: str -> Name Of Topology \\
        ------------------------------ \\
        RETURN: \\
        ---------------  \\
        toplist: list -> List Of Topology Files \\
        ------------------------------ \\
        '''

        toplist = []
        with open(topology_input) as ifile:
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^#include\s+\"(\S+)\"", line, flags=re.MULTILINE)
                if match:
                    match2 = re.search("ffbonded", match.group(1))
                    if match2:
                        continue
                    match2 = re.search("ffnonbonded", match.group(1))
                    if match2:
                        continue
                    match2 = re.search("forcefield.itp", match.group(1))
                    if match2:
                        continue
                    match2 = re.search("posre.itp", match.group(1))
                    if match2:
                        continue
                    foundname = match.group(1)
                    check = os.path.isfile(foundname)
                    if not check:
                        foundname = os.path.join(self.dict_input_userparameters['gmxtop_path'], *foundname.strip('/').split('/'))
                        check = os.path.isfile(foundname)
                        if not check:
                            print("File " + foundname + " was not found. Maybe update the gmxpath variable in the script? Exiting.")
                            exit(1)
                    toplist.append(foundname)

                    toplist.extend(self.get_list_topologies(foundname))
        return toplist

    def get_atoms(self, str_qmmm_topology) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        ---------------  \\
        Reads Elements For All Atoms \\
        ------------------------------ \\
        INPUT: \\
        ---------------  \\
        str_qmmm_topology: str -> Name Of QMMM Topology File \\
        ------------------------------ \\
        RETURN: \\
        ---------------  \\
        atoms: list -> List Of Element Of All Atoms \\
        ------------------------------ \\
        '''
        atoms = []
        str_file_mass_map = os.path.join('json_files', 'mass_map.json')
        with open(str_file_mass_map, 'r') as file:
            mass_map = json.load(file)

        name_map = {value: key for key, value in mass_map.items()}
        with open(str_qmmm_topology) as qmmm_topology_file:
            for line in qmmm_topology_file:
                match = re.search(r"\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    # logger(logfile, "moleculetype section was identified\n")
                    break
            for line in qmmm_topology_file:
                match = re.search(r"\[\s+atoms\s*\]", line, flags=re.MULTILINE)
                if match:
                    # logger(logfile, "atoms section was identified\n")
                    break
            for line in qmmm_topology_file:
                match = re.search(r"^\s*\[", line, flags=re.MULTILINE)
                if match:
                    break
                match = re.search(
                    r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-]*\d+[\.]*[\d+]*)\s+(\d+[\.]*[\d+]*)",
                    line,
                    flags=re.MULTILINE,
                )
                if match:
                    atomtype = str(match.group(2))
                    atommass = float(match.group(8))
                    foundname = ""
                    # find atom type based on mass
                    for key in name_map.items():
                        foundmass = key[0]
                        massdiff = math.sqrt(
                            (float(atommass) - float(foundmass))
                            * (float(atommass) - float(foundmass))
                        )
                        if massdiff < 0.05:
                            foundname = name_map[foundmass]
                            break
                    if foundname != "":
                        testmass = mass_map[foundname]
                        massdiff = math.sqrt(
                            (float(atommass) - float(foundmass))
                            * (float(atommass) - float(foundmass))
                        )
                        # XX change to logger
                        # if massdiff > 0.01:
                        #     logger(
                        #         logfile,
                        #         str(
                        #             "Found a mass of "
                        #             + str(atommass)
                        #             + " for atom type "
                        #             + str(atomtype)
                        #             + ' (identified as atom name "'
                        #             + str(foundname)
                        #             + '"), which is more than 0.01 different from the expected mass of '
                        #             + str(testmass)
                        #             + ". Atom index was "
                        #             + str(match.group(1))
                        #             + ". This has no effect except unless the atom was identified wrongly or dynamics are intended. Clean your ffnonbonded.itp to avoid these messages!\n"
                        #         ),
                        #     )
                        atoms.append(foundname)
                    else:
                        # logger(
                        #     logfile,
                        #     str(
                        #         "Atom type "
                        #         + str(atomtype)
                        #         + " could not be translated to a regular atom name. Exiting. Last line:\n"
                        #     ),
                        # )
                        # logger(logfile, line)
                        exit(1)
        return atoms
