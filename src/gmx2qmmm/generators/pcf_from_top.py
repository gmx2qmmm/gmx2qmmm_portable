# encoding: ISO-8859-15

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# This will take a gro and top file, and take the point charges from there.
# THIS IS NOT YET IMPLEMENTED#This list of charges and coordinates will then be trimmed by a set of point charges which then will be concatenated to a charge shift model, like in a QM/MM approach.
# THIS IS NOT YET IMPLEMENTED#The trimming is optional.
# will print output list as x y z charge in Angstrom
# THIS IS NOT YET TRUE#needs input pdb, file containing groups/atoms to trim (format: <number of groups>\n<for each group:><amount of atoms in group>\n<indices of atoms, one per line>
# needs input .gro file, input .top file, PCF output file

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$"  # during a rain storm

import re
import os
import sys


def checkformol(molname, inp):
    '''
    ------------------------------
    EFFECT: \\
    ---------------
    function checks if the molecule (molname) is listed in the topology file
    ------------------------------
    INPUT: \\
    ---------------
    molname: string, name of molecule
    inp: string, name of topology file
    ------------------------------
    RETURN: \\
    ---------------
    correct: -> bool
    ------------------------------
    '''
    with open(inp) as ifile:
        correct = False
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
                        matchstring = r"^\s*" + re.escape(molname)
                        match = re.search(matchstring, line, flags=re.MULTILINE)
                        if match:
                            correct = True
                        break
            if correct:
                break
    return correct

def getincludelist(inp, gmxtop_path):
    '''
    ------------------------------
    EFFECT: \\
    ---------------
    looks for other topology files referred to in the original file and adds them to the list of topology files
    ------------------------------
    INPUT: \\
    ---------------
    inp: string, name of topology file
    gmxtop_path: string, path to grolib
    ------------------------------
    RETURN: \\
    ---------------
    toplist: list of topology files
    ------------------------------
    '''
    toplist = []
    with open(inp) as ifile:
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
                    foundname = os.path.join(gmxtop_path, *foundname.strip('/').split('/'))
                    check = os.path.isfile(foundname)
                    if not check:
                        print("File " + foundname + " was not found. Maybe update the gmxpath variable in the script? Exiting.")
                        exit(1)
                toplist.append(foundname)

                toplist.extend(getincludelist(foundname, gmxtop_path))
    return toplist

def readcharges(molvecentry, top, gmxtop_path):
    '''
    ------------------------------
\\
    EFFECT: \\
    ---------------
\\
    reads charges for all atoms in a molecule type
\\
    ------------------------------
\\
    INPUT: \\
    ---------------
\\
    molvecentry: list with molecule name and amount of this molecule
    top: string, name of topology file
    gmxtop_path: string, path to grolib
    ------------------------------
\\
    RETURN: \\
    ---------------
\\
    finalvec: list of charges for every atom in this molecule types
\\
    ------------------------------
\\
    '''

    cvec = []
    curr_top = top
    molname = molvecentry[0]
    molcount = molvecentry[1]
    found = checkformol(molname, top)

    if not found:
        toplist = getincludelist(top, gmxtop_path)
        for element in toplist:
            found = checkformol(molname, element)
            if found:
                curr_top = element
                break
    if not found:
        print("No charges found for " + str(molname) + ". Exiting.")
        exit(1)
    with open(curr_top) as ifile:
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
                    matchstring = r"^\s*" + re.escape(molname)
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
                cvec.append(float(match.group(1)))

    finalcvec = []
    for i in range(0, int(molcount)):
        finalcvec.extend(cvec)
    return finalcvec

def readg96(inp):
    '''
    ------------------------------
    EFFECT: \\
    ---------------
    reads coordinates of all atoms from a g96 file
    ------------------------------
    INPUT: \\
    ---------------
    inp: string, name of g96 file
    ------------------------------
    RETURN: \\
    ---------------
    coords: list of coordinates
    ------------------------------
    '''
    coords = []
    with open(inp) as ifile:
        count = 0
        for line in ifile:
            count += 1
            if count == 4:
                break
        count = 1
        for line in ifile:
            match = re.search(
                r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                coords.append(float(match.group(5)) * 10.0)
                coords.append(float(match.group(6)) * 10.0)
                coords.append(float(match.group(7)) * 10.0)
            else:
                break
    return coords

def readgeo(inp):
    coords = []
    n_a = 0
    with open(inp) as ifile:
        for line in ifile:
            break
        for line in ifile:
            match = re.search(r"^\s*(\d+)", line, flags=re.MULTILINE)
            if match:
                n_a = int(match.group(1))
                break
            else:
                print(".gro is corrupt (no number of atoms found, second line). Exiting.")
                exit(1)
        count = 1
        for line in ifile:
            match = re.search(
                r"^(.{5})(.{5})(.{5})(.{5})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                coords.append(float(match.group(5)) * 10.0)
                coords.append(float(match.group(6)) * 10.0)
                coords.append(float(match.group(7)) * 10.0)
            else:
                print(".gro is corrupt. Exiting.")
                print("Last line:")
                print(line)
                exit(1)
            count += 1
            if count > n_a:
                break
    return coords

def readmols(top):
    '''
    ------------------------------
    EFFECT: \\
    ---------------
    reads list of molecules from the topology file
    ------------------------------
    INPUT: \\
    ---------------
    top: string, name of topology file
    ------------------------------
    RETURN: \\
    ---------------
    mollist: list of molecules and their amount in the system
    ------------------------------
    '''

    mollist = []
    with open(top) as ifile:
        found = False
        for line in ifile:
            match = re.search(r"^\[ molecules \]", line, flags=re.MULTILINE)
            if match:
                found = True
                break
        if not found:
            print('No "molecules" entry in ' + str(top) + " found. Exiting.")
            exit(1)
        for line in ifile:
            match = re.search(r"^;", line, flags=re.MULTILINE)
            if match:
                continue
            else:
                match = re.search(r"^(\S+)\s+(\d+)", line, flags=re.MULTILINE)
                if match:
                    mollist.append([match.group(1), match.group(2)])
                else:
                    print("Found an incomprehensible line in molecules list. Exiting.")
                    print("Last line was:")
                    print(line)
                    exit(1)
    return mollist

def read_numatoms(inp):
    '''
    ------------------------------
    EFFECT: \\
    ---------------
    reads and returns the number of atoms in the system
    ------------------------------
    INPUT: \\
    ---------------
    inp: string, name of structure file
    ------------------------------
    RETURN: \\
    ---------------
    int, number of atoms
    ------------------------------
    '''

    filetype = inp[-3:]
    if filetype == 'g96':
        pos_searched = False
        with open(inp ) as fp:
            for i, line in enumerate(fp):
                if 'POSITION' in line:
                    pos_start = i
                    pos_searched = True
                elif ('END' in line) and pos_searched:
                    pos_end = i
                    pos_searched = False
        return pos_end - (pos_start + 1)

    if filetype == 'gro':
        with open(inp) as fp:
            for i in range(2):
                num = fp.readline().split()

        return int(num[0])

    raise ValueError(f'Unknown file type {filetype}. Can only process ".gro" and ".g96"')


def makeout(coords, charges, name):
    with open(name, "w") as fp:
        for i in range(0, len(charges)):

            for j in range(0, 3):
                fp.write("{:>16.8f} ".format(coords[i * 3 + j]))
            fp.write("{:>16.8f}\n".format(charges[i]))


def generate_pcf_from_top(gro, top, out, gmxtop_path):
    chargevec = []
    mollist = readmols(top)
    for element in mollist:
        chargevec.extend(readcharges(element, top, gmxtop_path))
    term = str(str(gro[-3]) + str(gro[-2]) + str(gro[-1]))
    geo = []
    if term == "g96":
        geo = readg96(gro)
    else:
        geo = readgeo(gro)

    if len(geo) != 3 * len(chargevec):
        print("Not all atoms (" + str(
            len(geo) / 3.0
        ) + ") were replaced by charges (" + str(
            len(chargevec)
        ) + ") in full iteration step! Exiting.")
        exit(1)
    makeout(geo, chargevec, out)


if __name__ == "__main__":
    generate_pcf_from_top(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
