#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""
To change this license header, choose License Headers in Project
Properties. To change this template file, choose Tools | Templates and
open the template in the editor. THIS IS NOT ABLE TO RUN STANDALONE -
ONLY USE THE LISTSONLY FUNCTION This will take a gromacs .top file, and
remove all terms that are not needed for a QM/MM energy description.
Will result in a QM/MM viable .top file, if used with a proper energy
exlusion group set up. Will also set up an exclusion list between all QM
atoms. This is required in the MD run. Needs .top, file qmatomlist,
output file name, list of flags
"""

#   // MEATDATA //
__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$"  # During a rain storm

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os
import os.path
import numpy as np
from gmx2qmmm.generators import reindexing

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries

#   // TODOS & NOTES //
#   TODO:
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class GenerateTopology():

    '''
    This class generates a new topology file
    '''

    def __init__(self, input_dict) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Coordinates The Writing Of The New Large Topology File \\
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

        #self.system = system
        self.input_dict = input_dict
        #self.basedir = basedir
        self.mmflaglist = [] # figure out what that variable exactly is later (empty list in example) XX AJ
        self.qmmm_topology = str(input_dict['jobname'] + ".qmmm.top.init")
        self.generate_top_listsonly()

    def find_ffnonbonded(self, includedata):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Looks For ffnonbonded Files \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        ffnonbonded: str -> Name Of ffnonbonded File \\
        ------------------------------ \\
        '''

        ffnonbonded = ""
        for element in includedata:
            match = re.search("ffnonbonded\.itp", element, flags=re.MULTILINE)
            if match:
                if os.path.isfile(element):
                    ffnonbonded = element
                    break
            else:
                curr_dir = ""
                curr_dir_list = re.findall("\S+/", element)
                for entry in curr_dir_list:
                    curr_dir += entry
                new_ffnonbonded = str(curr_dir + "ffnonbonded.itp")
                if os.path.isfile(new_ffnonbonded):
                    ffnonbonded = new_ffnonbonded
                    break
        if ffnonbonded == "":
            pass
            # logger(
            #     logfile,
            #     str(
            #         "Did not find an ffnonbonded file. Check your masses in the qmmm.top file!\n"
            #     ),
            # )
        return ffnonbonded


    def get_mass(self, atomtype, ffnonbonded):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        XX AJ I'm not quite sure what this function is doing, it's not being called for my current example, check later \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        atomtype: str -> \\
        ffnonbonded: str -> Name Of ffnonbonded File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        minimal_mass: float -> XX? \\
        ------------------------------ \\
        '''

        matchstring = r"^\s*" + re.escape(atomtype) + "\s+\d+\s+(\d+\.\d+)"
        mass = []
        with open(ffnonbonded) as ifile:
            for line in ifile:
                match = re.search(matchstring, line, flags=re.MULTILINE)
                if match:
                    mass.append(float(match.group(1)))
        minimal_mass = min(mass)

        return minimal_mass

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
                        foundname = os.path.join(self.input_dict['gmxtop_path'], *foundname.strip('/').split('/'))
                        check = os.path.isfile(foundname)
                        if not check:
                            print("File " + foundname + " was not found. Maybe update the gmxpath variable in the script? Exiting.")
                            exit(1)
                    toplist.append(foundname)

                    toplist.extend(self.get_list_topologies(foundname))
        return toplist

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
        with open(self.input_dict['topologyfile'], 'r') as file_input:
            bool_match = False
            for line in file_input:
                match = re.search(r"^\[ molecules \]", line, flags=re.MULTILINE)
                if match:
                    bool_match = True
                    break
            if not bool_match:
                #   XX AJ turn into Logging
                print('No "molecules" entry in the provided topology found. Exiting.')
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

    def get_molecule_topology(self, molecule_name):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads The Topology Or ITP File For One Molecule \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        molecule_topologies: list -> List Of Sublists With Molecules And Their Topology Or ITP Files \\
        ------------------------------ \\
        '''

        foundtop = ""
        toplist = [self.input_dict['topologyfile']] + self.get_list_topologies(self.input_dict['topologyfile'])
        for element in toplist:
            with open(element) as ifile:
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
                            else:
                                matchstring = r"^\s*" + re.escape(molecule_name)
                                match = re.search(matchstring, line, flags=re.MULTILINE)
                                if match:
                                    foundtop = str(element)
                                    break
                            break
                    if foundtop != "":
                        break
                if foundtop != "":
                    break
            if foundtop != "":
                break
        if foundtop == "":
            print("Molecule " + str(molecule_name) + " was not found in any top. Exiting.")
            exit(1)
        return foundtop


    def get_all_molecule_topologies(self):
        # previous name "get_molfindlist"

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads The Topology Or ITP File For Every Molecule \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        molecule_topologies: list -> List Of Sublists With Molecules And Their Topology Or ITP Files \\
        ------------------------------ \\
        '''

        molecule_topologies = []

        for element in self.read_molecules():
            molecule_topology = self.get_molecule_topology(element[0])
            moltopline = [element[0], molecule_topology]
            molecule_topologies.append(moltopline)

        return molecule_topologies


    def get_mollength(self):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads The Amount Of Atoms For Each Molecule Type \\
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

        mollength = []
        for element in self.molfindlist:
            curr_length = 0
            with open(element[1]) as ifile:
                for line in ifile:
                    match = re.search(r"^\[\s*moleculetype\s*\]", line, flags=re.MULTILINE)
                    if match:
                        for line in ifile:
                            match = re.search(r"^;", line, flags=re.MULTILINE)
                            if match:
                                continue
                            else:
                                matchstring = r"^\s*" + re.escape(element[0])
                                match = re.search(matchstring, line, flags=re.MULTILINE)
                                if match:
                                    for line in ifile:
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
                                                match = re.search(
                                                    r"^\s*\n", line, flags=re.MULTILINE
                                                )
                                                if match:
                                                    break
                                                else:
                                                    curr_length += 1
                                            if curr_length > 0:
                                                break
                                    if curr_length > 0:
                                        break
                        if curr_length > 0:
                            break
            mollength.append(curr_length)
        return mollength


    def clean_exclusions(self, excludedata, number_of_atoms):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Cleans And Sorts The List For Exclucion Of Non-Bonded Interactions (XX AJ but I don't quite understand how) \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        excludedata: list -> List Of Atom Indices For Exclusion (XX?) \\
        number_of_atoms: int -> ? \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        final_excludedata: list -> List Of Atom Indices For Exclusion Sorted \\
        ------------------------------ \\
        '''

        # logger(logfile, str("Cleaning exclusion list...\n"))
        new_excludedata = []
        for i in range(1, number_of_atoms + 1):
            new_excludeline = [int(i)]
            new_excludedata.append(new_excludeline)
        for element in excludedata:
            for j in range(1, len(element)):
                if int(element[j]) > int(element[0]) and int(element[j]) not in np.array(
                    new_excludedata[int(element[0]) - 1]
                ).astype(int):
                    new_excludedata[int(element[0]) - 1].append(int(element[j]))
                if int(element[j]) < int(element[0]) and int(element[0]) not in np.array(
                    new_excludedata[int(element[j]) - 1]
                ).astype(int):
                    new_excludedata[int(element[j]) - 1].append(int(element[0]))
        final_excludedata = []
        for element in new_excludedata:
            if len(element) > 1:
                final_excludedata.append(element)
        # logger(logfile, str("Cleaning done.\n"))

        return final_excludedata


    def cleanagain_exclusions(self, excludedata): # used

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Cleans And Sorts The List For Exclucion Of Non-Bonded Interactions (XX AJ but I don't quite understand how, see clean_exclusions) \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        excludedata: list -> List Of Atom Indices For Exclusion (XX?) \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        new_excludedata: list -> List Of Atom Indices For Exclusion Sorted \\
        ------------------------------ \\
        '''

        # logger(logfile, str("Formatting exclusion list...\n"))
        new_excludedata = []
        for element in excludedata:
            for i in range(1, len(element)):
                new_excludedata.append([int(element[0]), int(element[i])])
        # logger(logfile, str("Formatting done.\n"))
        return new_excludedata

    def Permutate_Pairwise_no_rep(self, qmatomlist): 

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Creates a list of all unique pairwise permutations of the elements of a list
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        qmatomlist: list -> List Of qm atoms \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        pairwise_permutation: list -> List unique pairwise permuations of the qm atoms \\
        ------------------------------ \\
        '''

        qmatoms = qmatomlist
        pairwise_permutation = []
        for i in range(len(qmatoms)):
            for j in range(i+1, len(qmatoms)):
                pairwise_permutation.append([qmatoms[i], qmatoms[j]])

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

    def make_large_top(self):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Writes The Information From All Topology Files To One Large Topology \\
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

        mollength = self.get_mollength()
        red_molfindlist = []
        for entry in self.molfindlist:
            red_molfindlist.append(entry[1])
        offset = 0
        real_last_res = 0
        curr_res = 0
        real_last_cg = 0
        curr_cg = 0
        count = 0
        with open(self.qmmm_topology, "w") as ofile:
            includedata = []
            atomdata = []
            bonddata = []
            pairdata = []
            angledata = []
            dihedraldata = []
            settledata = []
            excludedata = []
            if self.input_dict['topologyfile'] not in red_molfindlist:
                with open(self.input_dict['topologyfile']) as ifile:
                    blocked = False
                    for line in ifile:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            ofile.write(line)
                            continue
                        match = re.search(r"^#ifdef\s+(\S+)", line, flags=re.MULTILINE)
                        if match and (match.group(1)).upper() not in self.mmflaglist:
                            blocked = True
                            continue
                        match = re.search(r"^#ifndef\s+(\S+)", line, flags=re.MULTILINE)
                        if match and (match.group(1)).upper() in self.mmflaglist:
                            blocked = True
                            continue
                        match = re.search(r"^#else", line, flags=re.MULTILINE)
                        if match and blocked:
                            blocked = False
                            continue
                        if match and not blocked:
                            blocked = True
                            continue
                        match = re.search(r"^#endif", line, flags=re.MULTILINE)
                        if match and blocked:
                            blocked = False
                        if blocked:
                            continue
                        match = re.search(
                            r"^#include\s+\"(\S+)\"", line, flags=re.MULTILINE
                        )
                        if match:
                            includedata.append(match.group(1))
                            continue
            for molecule in self.read_molecules():
                curr_topfile = ""
                for entry in self.molfindlist:
                    if molecule[0] == entry[0]:
                        curr_topfile = entry[1]
                        break
                for k in range(0, int(molecule[1])):
                    with open(curr_topfile) as ifile:
                        blocked = False
                        for line in ifile:
                            match = re.search(r"^;", line, flags=re.MULTILINE)
                            if match and k == 0:
                                ofile.write(line)
                                continue
                            match = re.search(r"^#ifdef\s+(\S+)", line, flags=re.MULTILINE)
                            if match and (match.group(1)).upper() not in self.mmflaglist:
                                blocked = True
                                continue
                            match = re.search(r"^#ifndef\s+(\S+)", line, flags=re.MULTILINE)
                            if match and (match.group(1)).upper() in self.mmflaglist:
                                blocked = True
                                continue
                            match = re.search(r"^#else", line, flags=re.MULTILINE)
                            if match and blocked:
                                blocked = False
                                continue
                            if match and not blocked:
                                blocked = True
                                continue
                            match = re.search(r"^#endif", line, flags=re.MULTILINE)
                            if match and blocked:
                                blocked = False
                            if blocked:
                                continue
                            match = re.search(
                                r"^#include\s+\"(\S+)\"", line, flags=re.MULTILINE
                            )
                            if match:
                                includedata.append(match.group(1))
                                continue
                            match = re.search(
                                r"^\[\s+moleculetype\s+\]", line, flags=re.MULTILINE
                            )
                            if match:
                                for line in ifile:
                                    matchstring = r"^\s*" + re.escape(molecule[0])
                                    match = re.search(matchstring, line, flags=re.MULTILINE)
                                    if match:
                                        blocked = False
                                        for line in ifile:
                                            match = re.search(
                                                r"^\[\s+moleculetype\s+\]",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if match:
                                                break
                                            match = re.search(
                                                r"^#ifdef\s+(\S+)", line, flags=re.MULTILINE
                                            )
                                            if (
                                                match
                                                and (match.group(1)).upper() not in self.mmflaglist
                                            ):
                                                blocked = True
                                                continue
                                            match = re.search(
                                                r"^#ifndef\s+(\S+)",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if (
                                                match
                                                and (match.group(1)).upper() in self.mmflaglist
                                            ):
                                                blocked = True
                                                continue
                                            match = re.search(
                                                r"^#else", line, flags=re.MULTILINE
                                            )
                                            if match and blocked:
                                                blocked = False
                                                continue
                                            if match and not blocked:
                                                blocked = True
                                                continue
                                            match = re.search(
                                                r"^#endif", line, flags=re.MULTILINE
                                            )
                                            if match and blocked:
                                                blocked = False
                                            if blocked:
                                                continue
                                            match = re.search(
                                                r"^#include\s+\"(\S+)\"",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if match:
                                                includedata.append(match.group(1))
                                                continue
                                            match = re.search(
                                                r"^\[\s+atoms\s+\]",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if match:
                                                for line in ifile:
                                                    # begin blockblock
                                                    match = re.search(
                                                        r"^#ifdef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        not in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#ifndef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#else", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                        continue
                                                    if match and not blocked:
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#endif", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                    if blocked:
                                                        continue
                                                    # end blockblock
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                                    match = re.search(
                                                        r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-]*\d+[\.]*[\d+]*)\s+(\d+[\.]*[\d]*)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if match:
                                                        if (
                                                            k == 0
                                                            and int(match.group(1)) == 1
                                                            or int(match.group(3))
                                                            != int(real_last_res)
                                                        ):
                                                            real_last_res = int(
                                                                match.group(3)
                                                            )
                                                            curr_res += 1
                                                        if (
                                                            k == 0
                                                            and int(match.group(1)) == 1
                                                            or int(match.group(6))
                                                            != int(real_last_cg)
                                                        ):
                                                            real_last_cg = int(
                                                                match.group(6)
                                                            )
                                                            curr_cg += 1
                                                        atomvec = [
                                                            int(match.group(1))
                                                            + int(offset),
                                                            match.group(2),
                                                            int(curr_res),
                                                            match.group(4),
                                                            match.group(5),
                                                            int(curr_cg),
                                                            match.group(7),
                                                            match.group(8),
                                                        ]
                                                        atomdata.append(atomvec)
                                                        continue
                                                    match = re.search(
                                                        r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-]*\d+[\.]*[\d]*)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if match:
                                                        if (
                                                            k == 0
                                                            and int(match.group(1)) == 1
                                                            or int(match.group(3))
                                                            != int(real_last_res)
                                                        ):
                                                            real_last_res = int(
                                                                match.group(3)
                                                            )
                                                            curr_res += 1
                                                        if (
                                                            k == 0
                                                            and int(match.group(1)) == 1
                                                            or int(match.group(6))
                                                            != int(real_last_cg)
                                                        ):
                                                            real_last_cg = int(
                                                                match.group(6)
                                                            )
                                                            curr_cg += 1
                                                        atomvec = [
                                                            int(match.group(1))
                                                            + int(offset),
                                                            match.group(2),
                                                            int(curr_res),
                                                            match.group(4),
                                                            match.group(5),
                                                            int(curr_cg),
                                                            match.group(7),
                                                        ]
                                                        atomdata.append(atomvec)
                                                real_last_cg = int(0)
                                                real_last_res = int(0)
                                                continue
                                            match = re.search(
                                                r"^\[\s*bonds\s*\]",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if match:
                                                for line in ifile:
                                                    # begin blockblock
                                                    match = re.search(
                                                        r"^#ifdef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        not in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#ifndef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#else", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                        continue
                                                    if match and not blocked:
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#endif", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                    if blocked:
                                                        continue
                                                    # end blockblock
                                                    match = re.search(
                                                        r"^;", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        continue
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                                    bondline = re.findall("\S+", line)
                                                    if bondline:
                                                        bondline[0] = int(
                                                            bondline[0]
                                                        ) + int(offset)
                                                        bondline[1] = int(
                                                            bondline[1]
                                                        ) + int(offset)
                                                        bonddata.append(bondline)
                                                continue
                                            match = re.search(
                                                r"^\[\s*pairs\s*\]",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if match:
                                                for line in ifile:
                                                    # begin blockblock
                                                    match = re.search(
                                                        r"^#ifdef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        not in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#ifndef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#else", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                        continue
                                                    if match and not blocked:
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#endif", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                    if blocked:
                                                        continue
                                                    # end blockblock
                                                    match = re.search(
                                                        r"^;", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        continue
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                                    pairline = re.findall("\S+", line)
                                                    if pairline:
                                                        pairline[0] = int(
                                                            pairline[0]
                                                        ) + int(offset)
                                                        pairline[1] = int(
                                                            pairline[1]
                                                        ) + int(offset)
                                                        pairdata.append(pairline)
                                                continue
                                            match = re.search(
                                                r"^\[\s*angles\s*\]",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if match:
                                                for line in ifile:
                                                    # begin blockblock
                                                    match = re.search(
                                                        r"^#ifdef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        not in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#ifndef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#else", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                        continue
                                                    if match and not blocked:
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#endif", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                    if blocked:
                                                        continue
                                                    # end blockblock
                                                    match = re.search(
                                                        r"^;", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        continue
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                                    angleline = re.findall("\S+", line)
                                                    if angleline:
                                                        angleline[0] = int(
                                                            angleline[0]
                                                        ) + int(offset)
                                                        angleline[1] = int(
                                                            angleline[1]
                                                        ) + int(offset)
                                                        angleline[2] = int(
                                                            angleline[2]
                                                        ) + int(offset)
                                                        angledata.append(angleline)
                                                continue
                                            match = re.search(
                                                r"^\[ dihedrals \]",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if match:
                                                for line in ifile:
                                                    # begin blockblock
                                                    match = re.search(
                                                        r"^#ifdef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        not in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#ifndef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#else", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                        continue
                                                    if match and not blocked:
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#endif", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                    if blocked:
                                                        continue
                                                    # end blockblock
                                                    match = re.search(
                                                        r"^;", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        continue
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                                    dihedralline = re.findall("\S+", line)
                                                    if dihedralline:
                                                        dihedralline[0] = int(
                                                            dihedralline[0]
                                                        ) + int(offset)
                                                        dihedralline[1] = int(
                                                            dihedralline[1]
                                                        ) + int(offset)
                                                        dihedralline[2] = int(
                                                            dihedralline[2]
                                                        ) + int(offset)
                                                        dihedralline[3] = int(
                                                            dihedralline[3]
                                                        ) + int(offset)
                                                        dihedraldata.append(dihedralline)
                                                continue
                                            match = re.search(
                                                r"^\[ settles \]", line, flags=re.MULTILINE
                                            )
                                            if match:
                                                for line in ifile:
                                                    # begin blockblock
                                                    match = re.search(
                                                        r"^#ifdef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        not in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#ifndef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#else", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                        continue
                                                    if match and not blocked:
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#endif", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                    if blocked:
                                                        continue
                                                    # end blockblock
                                                    match = re.search(
                                                        r"^;", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        continue
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                                    settleline = re.findall("\S+", line)
                                                    if settleline:
                                                        settleline[0] = int(
                                                            settleline[0]
                                                        ) + int(offset)
                                                        settledata.append(settleline)
                                                continue
                                            match = re.search(
                                                r"^\[ exclusions \]",
                                                line,
                                                flags=re.MULTILINE,
                                            )
                                            if match:
                                                for line in ifile:
                                                    # begin blockblock
                                                    match = re.search(
                                                        r"^#ifdef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        not in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#ifndef\s+(\S+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if (
                                                        match
                                                        and (match.group(1)).upper()
                                                        in self.mmflaglist
                                                    ):
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#else", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                        continue
                                                    if match and not blocked:
                                                        blocked = True
                                                        continue
                                                    match = re.search(
                                                        r"^#endif", line, flags=re.MULTILINE
                                                    )
                                                    if match and blocked:
                                                        blocked = False
                                                    if blocked:
                                                        continue
                                                    # end blockblock
                                                    match = re.search(
                                                        r"^;", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        continue
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                                    excludeline = re.findall("\S+", line)
                                                    if excludeline:
                                                        for i in range(0, len(excludeline)):
                                                            excludeline[i] = int(
                                                                excludeline[i]
                                                            ) + int(offset)
                                                        excludedata.append(excludeline)
                                                continue
                    offset += mollength[count]
                count += 1
            for i in range(0, len(includedata)):
                ofile.write('#include "' + str(includedata[i]) + '"\n')
            ofile.write("\n[ moleculetype ]\nQMMM_model     3\n\n[ atoms ]\n")
            ffnb = self.find_ffnonbonded(includedata)

            if self.input_dict['useinnerouter']:
                outeratomlist = set([int(i) for i in self.read_atoms_list(self.input_dict['inneratomslist'])])
                atomdatanew = []
                bonddatanew = []
                pairdatanew = []
                angledatanew = []
                dihedraldatanew = []
                settledatanew = []
                excludedatanew = []
                mem = [] #[[old index, new index]]
                c = 0
                #Create a memory container and update atomdata
                for j in range(len(atomdata) - 1):
                    if int(atomdata[j][0]) in outeratomlist:
                        mem.append([atomdata[j + 1][0], atomdata[j][0]])
                        if mem[-1][0] != mem[-1][1]:
                            c += 1
                            mem[-1] = [mem[-1][1], mem[-1][0] - c]
                        else:
                            continue
                    else:
                        atomdatanew.append(atomdata[j])
                        if not mem == [] and mem[-1][0] != mem[-1][1]:
                            mem.append([atomdata[j + 1][0], atomdata[j][0]])
                            if mem[-1][0] != mem[-1][1]:
                                mem[-1] = [mem[-1][1], mem[-1][0] - (c + 1)]
        
                            atomdatanew[-1][0] = mem[-1][1]
                            atomdatanew[-1][5] = mem[-1][1]
                        else:
                            mem.append([atomdata[j][0], atomdata[j][0]])
                            continue
                
                if int(atomdata[-1][0]) not in outeratomlist:
                    if int(atomdata[-2][0]) not in outeratomlist:
                        mem.append([int(atomdata[-1][0]), mem[-1][1] + 1])
                    else:
                        mem.append([int(atomdata[-1][0]), mem[-1][1]])

                    atomdatanew.append(atomdata[-1])
                    atomdatanew[-1][0] = mem[-1][1]
                    atomdatanew[-1][5] = mem[-1][1]
                else:
                    mem.append([int(atomdata[-1][0]), mem[-1][1]])

                '''if atomdatanew[-1][0] == atomdatanew[-2][0] and atomdatanew[-1][5] == atomdatanew[-2][5]:
                    atomdatanew[-1][0] += 1
                    atomdatanew[-1][5] += 1  
                    mem[-1][1] += 1  
                for i in range(1, len(atomdatanew)):
                    if atomdatanew[i][0] <= atomdatanew[i - 1][0]:
                        atomdatanew[i][0] = atomdatanew[i - 1][0] + 1
                        atomdatanew[i][5] = atomdatanew[i][0]
                for i in range(1, len(mem)):
                    if mem[i][1] <= mem[i - 1][1]:
                        mem[i][1] = mem[i - 1][1] + 1  ''' 

                c = 1
                for i in range(len(atomdatanew)):
                    atomdatanew[i][0] = c
                    atomdatanew[i][5] = c
                    c += 1
                memory = [] #[[old index, new index]]
                for i in mem:
                    if int(i[0]) not in outeratomlist:
                        memory.append(i)
                    else:
                        continue
            
                c = 1
                for i in range(len(memory)):
                    memory[i][1] = c
                    c += 1
                
                memory_dict = {item[0]: item[1] for item in memory}

                # Ensure that residue indecies in atomdatanew are incremental

                residue_map = {}
                curr = 1

                for i in range(len(atomdatanew)):
                    if int(atomdatanew[i][2]) not in residue_map:
                        residue_map[int(atomdatanew[i][2])] = int(curr)
                        curr += 1

                for i in range(len(atomdatanew)):
                    atomdatanew[i][2] = residue_map[int(atomdatanew[i][2])]

                #Update bonddata
                
                for j in range(len(bonddata)):
                    if (int(bonddata[j][0]) in outeratomlist) or (int(bonddata[j][1]) in outeratomlist):
                        continue
                    else:
                        bonddatanew.append(bonddata[j])
                        bonddatanew[-1][0] = memory_dict[bonddatanew[-1][0]]
                        bonddatanew[-1][1] = memory_dict[bonddatanew[-1][1]]

                #Update pairdata
                for j in range(len(pairdata)):
                    if (int(pairdata[j][0]) in outeratomlist) or (int(pairdata[j][1]) in outeratomlist):
                        continue
                    else:
                        pairdatanew.append(pairdata[j])
                        pairdatanew[-1][0] = memory_dict[pairdatanew[-1][0]]
                        pairdatanew[-1][1] = memory_dict[pairdatanew[-1][1]]
                
                #Update angledata
                
                for j in range(len(angledata)):
                    if (int(angledata[j][0]) in outeratomlist) or (int(angledata[j][1]) in outeratomlist) or (int(angledata[j][2]) in outeratomlist):
                        continue
                    else:
                        angledatanew.append(angledata[j])
                        angledatanew[-1][0] = memory_dict[angledatanew[-1][0]]
                        angledatanew[-1][1] = memory_dict[angledatanew[-1][1]]
                        angledatanew[-1][2] = memory_dict[angledatanew[-1][2]]

                #Update dihedraldata

                for j in range(len(dihedraldata)):
                    if (int(dihedraldata[j][0]) in outeratomlist) or (int(dihedraldata[j][1]) in outeratomlist) or (int(dihedraldata[j][2]) in outeratomlist) or (int(dihedraldata[j][3]) in outeratomlist):
                        continue
                    else:
                        dihedraldatanew.append(dihedraldata[j])
                        dihedraldatanew[-1][0] = memory_dict[dihedraldatanew[-1][0]]
                        dihedraldatanew[-1][1] = memory_dict[dihedraldatanew[-1][1]]
                        dihedraldatanew[-1][2] = memory_dict[dihedraldatanew[-1][2]]
                        dihedraldatanew[-1][3] = memory_dict[dihedraldatanew[-1][3]]

                #Update settledata

                for j in range(len(settledata)):
                    if int(settledata[j][0]) in outeratomlist:
                        continue
                    else:
                        settledatanew.append(settledata[j])
                        settledatanew[-1][0] = memory_dict[settledatanew[-1][0]]

                #Update excludedata

                #for j in range(len(excludedata)):
                #    if (int(excludedata[j][0]) in outeratomlist) or (int(excludedata[j][1]) in outeratomlist) or (int(excludedata[j][2]) in outeratomlist):
                #        continue
                #    else:
                #        excludedatanew.append(excludedata[j])
                #        excludedatanew[-1][0] = memory_dict[excludedatanew[-1][0]]
                #        excludedatanew[-1][1] = memory_dict[excludedatanew[-1][1]]
                for excl in excludedata:
                    if any(atom in outeratomlist for atom in excl):
                        continue
                    else:
                        new_exclusion = [memory_dict[atom] for atom in excl]
                        excludedatanew.append(new_exclusion)
   
                atomdata = atomdatanew
                bonddata = bonddatanew
                pairdata = pairdatanew
                angledata = angledatanew
                dihedraldata = dihedraldatanew
                settledata = settledatanew
                excludedata = excludedatanew

                #Re-index the qmatomlist based on memory_dict - Unnesessary, done in make_exclude_list globally
                #for i in range(len(qmatomlist)):
                #    qmatomlist[i] = memory_dict[qmatomlist[i]]

                qmatomlist = set([int(i) for i in self.read_atoms_list(self.input_dict['qmatomslist'])])
                for element in atomdata:
                    if int(element[0]) in qmatomlist:
                        ofile.write(
                            "{:>6d} {:>10s} {:>6d} {:>6s} {:>6s} {:>6d} {:>10.4f}".format(
                                int(element[0]),
                                str(element[1]),
                                int(element[2]),
                                str(element[3]),
                                str(element[4]),
                                int(element[5]),
                                float(0.0),
                            )
                        )
                    else:
                        ofile.write(
                            "{:>6d} {:>10s} {:>6d} {:>6s} {:>6s} {:>6d} {:>10.4f}".format(
                                int(element[0]),
                                str(element[1]),
                                int(element[2]),
                                str(element[3]),
                                str(element[4]),
                                int(element[5]),
                                float(element[6]),
                            )
                        )
                    if len(element) > 7:
                        ofile.write(" {:>10s}".format(str(element[7])))
                    else:
                        ofile.write(
                            " {:>10s}".format(str(self.get_mass(str(element[1]), ffnb)))
                        )
                    ofile.write("\n")
                ofile.write("\n[ bonds ]\n")

                for element in bonddata:
                    if (int(element[0]) in qmatomlist) or (
                        int(element[1]) in qmatomlist
                    ):
                        excludeline = [element[0], element[1]]
                        excludedata.append(excludeline)
                        continue
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
            
            
                '''for element in angledata:
                    if (
                        (int(element[0]) in qmatomlist)
                        and (int(element[1]) in qmatomlist)
                    ) or (
                        (int(element[1]) in qmatomlist)
                        and (int(element[2]) in qmatomlist)
                    ):
                        excludeline = [element[0], element[2]]
                        excludedata.append(excludeline)
                    else:
                        continue
                
                for element in dihedraldata:
                    if (
                        (int(element[0]) in qmatomlist)
                        and (int(element[1]) in qmatomlist)
                        and (int(element[2]) in qmatomlist)
                    ) or (
                        (int(element[1]) in qmatomlist)
                        and (int(element[2]) in qmatomlist)
                        and (int(element[3]) in qmatomlist)
                    ):
                        excludeline = [element[0], element[3]]
                        excludedata.append(excludeline)
                    else:
                        continue
                qmatomlist = list(qmatomlist)
                qmatomlist.sort()
                # increase exclusions for each Q-M1 atom
                for element in self.system.list_atoms_m1:
                    for entry in qmatomlist:
                        excludedata.append([int(element), int(entry)])
                # add other link correction exclusions: m2-q1,2, m3-q1
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        excludedata.append([int(self.system.list_atoms_m2[i][j]), int(self.system.list_atoms_q1[i])])
                        for k in range(0, len(self.system.list_atoms_q2[i])):
                            excludedata.append([int(self.system.list_atoms_m2[i][j]), int(self.system.list_atoms_q2[i][k])])
                for i in range(0, len(self.system.list_atoms_m3)):
                    for j in range(0, len(self.system.list_atoms_m3[i])):
                        for k in range(0, len(self.system.list_atoms_m3[i][j])):
                            excludedata.append([int(self.system.list_atoms_m3[i][j][k]), int(self.system.list_atoms_q1[i])])
                excludedata = self.clean_exclusions(excludedata, offset)
                excludedata = self.cleanagain_exclusions(excludedata)
                for i in Permutate_Pairwise_no_rep(qmatomlist): #Added by Nicola
                    if i not in excludedata and i.reverse() not in excludedata: #Added by Nicola
                        excludedata.append(i) #Added by Nicola
                    else: #Added by Nicola
                        continue #Added by Nicola
                qmatomlist = set([int(i) for i in qmatomlist])'''

                ofile.write("\n[ pairs ]\n")
                for element in pairdata:
                    if ([element[0], element[1]] or [element[1], element[0]]) not in excludedata: #Added by Nicola
                        for entry in element:
                            ofile.write(str(entry) + " ")
                    else:
                        continue
                    ofile.write("\n")
                ofile.write("\n[ angles ]\n")
                for element in angledata:
                    if (
                        (int(element[0]) in qmatomlist)
                        and (int(element[1]) in qmatomlist)
                    ) or (
                        (int(element[1]) in qmatomlist)
                        and (int(element[2]) in qmatomlist)
                    ):
                        #excludeline = [element[0], element[2]]
                        #excludedata.append(excludeline)
                        continue
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
                ofile.write("\n[ dihedrals ]\n")
                for element in dihedraldata:
                    if (
                        (int(element[0]) in qmatomlist)
                        and (int(element[1]) in qmatomlist)
                        and (int(element[2]) in qmatomlist)
                    ) or (
                        (int(element[1]) in qmatomlist)
                        and (int(element[2]) in qmatomlist)
                        and (int(element[3]) in qmatomlist)
                    ):
                        #excludeline = [element[0], element[3]]
                        #excludedata.append(excludeline)
                        continue
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
                ofile.write("\n[ settles ]\n")
                for element in settledata:
                    if int(element[0]) in qmatomlist:
                        continue
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
            
            
                ofile.write("\n[ exclusions ]\n")
                for element in excludedata:
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
                ofile.write("\n[ system ]\nProtein\n\n[ molecules ]\nQMMM_model 1")

            else:
                for element in atomdata:
                    if int(element[0]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int):
                        ofile.write(
                            "{:>6d} {:>10s} {:>6d} {:>6s} {:>6s} {:>6d} {:>10.4f}".format(
                                int(element[0]),
                                str(element[1]),
                                int(element[2]),
                                str(element[3]),
                                str(element[4]),
                                int(element[5]),
                                float(0.0),
                            )
                        )
                    else:
                        ofile.write(
                            "{:>6d} {:>10s} {:>6d} {:>6s} {:>6s} {:>6d} {:>10.4f}".format(
                                int(element[0]),
                                str(element[1]),
                                int(element[2]),
                                str(element[3]),
                                str(element[4]),
                                int(element[5]),
                                float(element[6]),
                            )
                        )
                    if len(element) > 7:
                        ofile.write(" {:>10s}".format(str(element[7])))
                    else:
                        ofile.write(
                            " {:>10s}".format(str(self.get_mass(str(element[1]), ffnb)))
                        )
                    ofile.write("\n")
                ofile.write("\n[ bonds ]\n")
                for element in bonddata:
                    if (int(element[0]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int)) or (
                        int(element[1]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int)
                    ):
                        excludeline = [element[0], element[1]]
                        excludedata.append(excludeline)
                        continue
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
                ofile.write("\n[ pairs ]\n")
                for element in pairdata:
                    if ([element[0], element[1]] or [element[1], element[0]]) not in excludedata: #Added by Nicola
                        for entry in element:
                            ofile.write(str(entry) + " ")
                    else:
                        continue
                    ofile.write("\n")
                ofile.write("\n[ angles ]\n")
                for element in angledata:
                    if (
                        (int(element[0]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                        and (int(element[1]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                    ) or (
                        (int(element[1]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                        and (int(element[2]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                    ):
                        excludeline = [element[0], element[2]]
                        excludedata.append(excludeline)
                        continue
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
                ofile.write("\n[ dihedrals ]\n")
                for element in dihedraldata:
                    if (
                        (int(element[0]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                        and (int(element[1]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                        and (int(element[2]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                    ) or (
                        (int(element[1]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                        and (int(element[2]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                        and (int(element[3]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int))
                    ):
                        excludeline = [element[0], element[3]]
                        excludedata.append(excludeline)
                        continue
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
                ofile.write("\n[ settles ]\n")
                for element in settledata:
                    if int(element[0]) in np.array(self.read_atoms_list(self.input_dict['qmatomslist'])).astype(int):
                        continue
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
                '''# increase exclusions for each Q-M1 atom
                for element in self.system.list_atoms_m1:
                    for entry in self.system.list_atoms_qm:
                        excludedata.append([int(element), int(entry)])
                # add other link correction exclusions: m2-q1,2, m3-q1
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        excludedata.append([int(self.system.list_atoms_m2[i][j]), int(self.system.list_atoms_q1[i])])
                        for k in range(0, len(self.system.list_atoms_q2[i])):
                            excludedata.append([int(self.system.list_atoms_m2[i][j]), int(self.system.list_atoms_q2[i][k])])
                for i in range(0, len(self.system.list_atoms_m3)):
                    for j in range(0, len(self.system.list_atoms_m3[i])):
                        for k in range(0, len(self.system.list_atoms_m3[i][j])):
                            excludedata.append([int(self.system.list_atoms_m3[i][j][k]), int(self.system.list_atoms_q1[i])])

            # XX AJ maybe combine these two functions? I don't quite understand how they work though
                excludedata = self.clean_exclusions(excludedata, offset)
                excludedata = self.cleanagain_exclusions(excludedata)
                for i in self.Permutate_Pairwise_no_rep(self.system.list_atoms_qm): 
                    if i not in excludedata and i.reverse() not in excludedata: 
                        excludedata.append(i)
                    else:
                        continue'''

                ofile.write("\n[ exclusions ]\n")
                for element in excludedata:
                    for entry in element:
                        ofile.write(str(entry) + " ")
                    ofile.write("\n")
                ofile.write("\n[ system ]\nProtein\n\n[ molecules ]\nQMMM_model 1")

        return None


    '''def make_gmx_index_file(self):
        # previous name "make_exclude_index" XX

        ''''''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Creates An Index File For Gromacs With The Different Atom Groups \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        ''''''

        outname = str(self.qmmm_topology) + ".ndx"
        with open(outname, "w") as ofile:
            count = 0
            ofile.write("[ QM ]\n")
            count = 0
            for element in self.system.list_atoms_qm:
                count += 1
                ofile.write(str(int(element)) + " ")
                if count == 15:
                    ofile.write("\n")
                    count = 0
            ofile.write("\n")
            if not self.system.list_atoms_inner == []:
                ofile.write("[ INNER ]\n")
                count=0
                for element in self.system.list_atoms_inner:
                    count += 1
                    ofile.write(str(int(element)) + " ")
                    if count == 15:
                        ofile.write("\n")
                        count = 0
                ofile.write("\n")
            if not self.system.list_atoms_outer == []:
                ofile.write("[ OUTER ]\n")
                count=0
                for element in self.system.list_atoms_outer:
                        count+=1
                        ofile.write(str(int(element)) + " ")
                        if count==15:
                                ofile.write("\n")
                                count=0
                ofile.write("\n")

        return None'''


    def generate_top_listsonly(self):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Coordinates Functions To Write The Index File For Gromacs And The Topology File \\
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

        self.molfindlist = self.get_all_molecule_topologies()
        #self.make_gmx_index_file()
        self.make_large_top()

        return None


    # def generate_top(self, sysargs):
    #     # XX AJ this function is only called when file is called as main, adapt this function later
    #     basedir = os.path.dirname(os.path.abspath(__file__))
    #     top = sysargs[1]
    #     qmatoms = sysargs[2]
    #     outname = sysargs[3]
    #     flaglist = sysargs[4]
    #     m1list = sysargs[5]
    #     qmatomlist = prep_pcf.read_qmatom_list(qmatoms)
    #     mollist = make_pcf.readmols(top)
    #     includelist = make_pcf.getincludelist(top, pathinfo)
    #     molfindlist = get_molfindlist(mollist, top, includelist)
    #     make_large_top()




    if __name__ == "__main__":
        import sys

        # generate_top(sys.argv)
