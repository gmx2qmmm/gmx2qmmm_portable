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

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$"  # during a rain storm

import os
import re
import numpy as np

from gmx2qmmm.pointcharges import generate_pcf_from_top as make_pcf
from gmx2qmmm.pointcharges import prepare_pcf_for_shift as prep_pcf
from gmx2qmmm.operations import reindexing # Added by Nicola


def logger(log, logstring):
    from datetime import datetime

    with open(log, "a") as ofile:
        ofile.write(str(datetime.now()) + " " + logstring)


def find_ffnonbonded(includedata, logfile):
    import os.path
    import re

    ffnb = ""
    for element in includedata:
        match = re.search("ffnonbonded\.itp", element, flags=re.MULTILINE)
        if match:
            if os.path.isfile(element):
                ffnb = element
                break
        else:
            curr_dir = ""
            curr_dir_list = re.findall("\S+/", element)
            for entry in curr_dir_list:
                curr_dir += entry
            new_ffnb = str(curr_dir + "ffnonbonded.itp")
            if os.path.isfile(new_ffnb):
                ffnb = new_ffnb
                break
    if ffnb == "":
        logger(
            logfile,
            str(
                "Did not find an ffnonbonded file. Check your masses in the qmmm.top file!\n"
            ),
        )
    return ffnb


def get_mass(atomtype, ffnb, logfile):
    import re

    matchstring = r"^\s*" + re.escape(atomtype) + "\s+\d+\s+(\d+\.\d+)"
    mass = []
    with open(ffnb) as ifile:
        for line in ifile:
            match = re.search(matchstring, line, flags=re.MULTILINE)
            if match:
                mass.append(float(match.group(1)))
    return min(mass)


def molfinder(top, includelist, molname):
    import re

    foundtop = ""
    toplist = [top] + includelist
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
                            matchstring = r"^\s*" + re.escape(molname)
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
        print("Molecule " + str(molname) + " was not found in any top. Exiting.")
        exit(1)
    return foundtop


def get_molfindlist(mollist, top, includelist):
    molfindlist = []
    for element in mollist:
        topfound = molfinder(top, includelist, element[0])
        moltopline = [element[0], topfound]
        molfindlist.append(moltopline)
    return molfindlist


def get_full_include_list(top, includelist):
    import re

    toplist = [top] + includelist
    full_include_list = [str(top)]
    for element in toplist:
        with open(element) as ifile:
            for line in ifile:
                match = re.search(r"^#include\s+\"(\S+)\"", line, flags=re.MULTILINE)
                if match:
                    if str(match.group(1)) not in full_include_list and (
                        str(match.group(1)) == "posre.itp"
                    ):
                        full_include_list.append(str(match.group(1)))
    return full_include_list


def get_mollength(molfindlist):
    import re

    mollength = []
    for element in molfindlist:
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


def clean_exclusions(excludedata, n_a, logfile):
    logger(logfile, str("Cleaning exclusion list...\n"))
    new_excludedata = []
    for i in range(1, n_a + 1):
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
    logger(logfile, str("Cleaning done.\n"))
    return final_excludedata


def cleanagain_exclusions(excludedata, logfile):
    logger(logfile, str("Formatting exclusion list...\n"))
    new_excludedata = []
    for element in excludedata:
        for i in range(1, len(element)):
            new_excludedata.append([int(element[0]), int(element[i])])
    logger(logfile, str("Formatting done.\n"))
    return new_excludedata

def Permutate_Pairwise_no_rep(qmatomlist): #Added by Nicola

    qmatoms = qmatomlist
    pairwise_permutation = []
    for i in range(len(qmatoms)):
        for j in range(i+1, len(qmatoms)):
            pairwise_permutation.append([qmatoms[i], qmatoms[j]])

    return pairwise_permutation

def Permutate_Pairwise(inneratomlist, outeratomlist): #Added by Nicola

    inneratoms = inneratomlist
    outeratoms = outeratomlist
    pairwise_permutation = []
    for i in range(len(inneratoms)):
        for j in range(len(outeratoms)):
            pairwise_permutation.append([inneratoms[i], outeratoms[j]])

    return pairwise_permutation

def merge_sublists_by_first_element(data): #Added by Nicola
    grouped = {}
    
    # Group sublists by their starting number
    for sublist in data:
        if sublist:  # Ensure the sublist is not empty
            key = sublist[0]
            if key not in grouped:
                grouped[key] = sublist[1:]  # Create a new entry
            else:
                grouped[key].extend(sublist[1:])  # Extend the existing entry

    # Create the merged list of sublists
    merged = [[key] + values for key, values in grouped.items()]
    
    return merged

def get_keys_from_value(dictionary, target_value): #Added by Nicola
    keys = []
    for key, value in dictionary.items():
        if value == target_value:
            keys.append(key)
    return keys




def make_large_top(
    top,
    molfindlist,
    mollist,
    qmatomlist,
    includelist,
    outname,
    flaglist,
    logfile,
    inneratomlist=[],
    outeratomlist=[],
    gmxplus = True 
): #gmxplus = True inneratomslist outeratomlist Added by Nicola
    mollength = get_mollength(molfindlist)
    red_molfindlist = []
    for entry in molfindlist:
        red_molfindlist.append(entry[1])
    offset = 0
    real_last_res = 0
    curr_res = 0
    real_last_cg = 0
    curr_cg = 0
    count = 0
    with open(outname, "w") as ofile:
        includedata = []
        atomdata = []
        bonddata = []
        pairdata = []
        angledata = []
        dihedraldata = []
        settledata = []
        excludedata = []
        if top not in red_molfindlist:
            with open(top) as ifile:
                blocked = False
                for line in ifile:
                    match = re.search(r"^;", line, flags=re.MULTILINE)
                    if match:
                        ofile.write(line)
                        continue
                    match = re.search(r"^#ifdef\s+(\S+)", line, flags=re.MULTILINE)
                    if match and (match.group(1)).upper() not in flaglist:
                        blocked = True
                        continue
                    match = re.search(r"^#ifndef\s+(\S+)", line, flags=re.MULTILINE)
                    if match and (match.group(1)).upper() in flaglist:
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
        for molecule in mollist:
            curr_topfile = ""
            for entry in molfindlist:
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
                        if match and (match.group(1)).upper() not in flaglist:
                            blocked = True
                            continue
                        match = re.search(r"^#ifndef\s+(\S+)", line, flags=re.MULTILINE)
                        if match and (match.group(1)).upper() in flaglist:
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
                                            and (match.group(1)).upper() not in flaglist
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
                                            and (match.group(1)).upper() in flaglist
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
                                                    not in flaglist
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
                                                    in flaglist
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
                                                    not in flaglist
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
                                                    in flaglist
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
                                                    not in flaglist
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
                                                    in flaglist
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
                                                    not in flaglist
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
                                                    in flaglist
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
                                                    not in flaglist
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
                                                    in flaglist
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
                                                    not in flaglist
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
                                                    in flaglist
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
                                                    not in flaglist
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
                                                    in flaglist
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
        ffnb = find_ffnonbonded(includedata, logfile)
        
        if not outeratomlist == []:        
            outeratomlist = set([int(i) for i in outeratomlist])
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
            for j in range(len(atomdata) -1 ):
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
                    mem[i][1] = mem[i - 1][1] + 1

            used_indices = set()
            last_index = 0
            for atom in atomdatanew:
                if atom[0] in used_indices or atom[0] <= last_index:
                    atom[0] = last_index + 1
                    atom[5] = atom[0]
                used_indices.add(atom[0])
                last_index = atom[0]'''
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

            for j in range(len(excludedata)):
                if (int(excludedata[j][0]) in outeratomlist) or (int(excludedata[j][1]) in outeratomlist):
                    continue
                else:
                    excludedatanew.append(excludedata[j])
                    excludedatanew[-1][0] = memory_dict[excludedatanew[-1][0]]
                    excludedatanew[-1][1] = memory_dict[excludedatanew[-1][1]]
                
            atomdata = atomdatanew
            bonddata = bonddatanew
            pairdata = pairdatanew
            angledata = angledatanew
            dihedraldata = dihedraldatanew
            settledata = settledatanew
            excludedata = excludedatanew

            #Re-index the qmatomlist based on memory_dict - Unnesessary, done in make_exclude_list globally
            for i in range(len(qmatomlist)):
                qmatomlist[i] = memory_dict[qmatomlist[i]]


            qmatomlist = set([int(i) for i in qmatomlist])
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
                        " {:>10s}".format(str(get_mass(str(element[1]), ffnb, logfile)))
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
            ofile.write("\n[ pairs ]\n")
            for element in pairdata:
                for entry in element:
                    ofile.write(str(entry) + " ")
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
                    excludeline = [element[0], element[2]]
                    excludedata.append(excludeline)
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
                    excludeline = [element[0], element[3]]
                    excludedata.append(excludeline)
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
                if int(element[0]) in np.array(qmatomlist).astype(int):
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
                        " {:>10s}".format(str(get_mass(str(element[1]), ffnb, logfile)))
                    )
                ofile.write("\n")
            ofile.write("\n[ bonds ]\n")

            for element in bonddata:
                if (int(element[0]) in np.array(qmatomlist).astype(int)) or (
                    int(element[1]) in np.array(qmatomlist).astype(int)
                ):
                    excludeline = [element[0], element[1]]
                    excludedata.append(excludeline)
                    continue
                for entry in element:
                    ofile.write(str(entry) + " ")
                ofile.write("\n")
            ofile.write("\n[ pairs ]\n")
            for element in pairdata:
                for entry in element:
                    ofile.write(str(entry) + " ")
                ofile.write("\n")
            ofile.write("\n[ angles ]\n")
            for element in angledata:
                if (
                    (int(element[0]) in np.array(qmatomlist).astype(int))
                    and (int(element[1]) in np.array(qmatomlist).astype(int))
                ) or (
                    (int(element[1]) in np.array(qmatomlist).astype(int))
                    and (int(element[2]) in np.array(qmatomlist).astype(int))
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
                    (int(element[0]) in np.array(qmatomlist).astype(int))
                    and (int(element[1]) in np.array(qmatomlist).astype(int))
                    and (int(element[2]) in np.array(qmatomlist).astype(int))
                ) or (
                    (int(element[1]) in np.array(qmatomlist).astype(int))
                    and (int(element[2]) in np.array(qmatomlist).astype(int))
                    and (int(element[3]) in np.array(qmatomlist).astype(int))
                ):
                    excludeline = [element[0], element[3]]
                    excludedata.append(excludeline)
                    continue
                for entry in element:
                    ofile.write(str(entry) + " ")
                ofile.write("\n")
            ofile.write("\n[ settles ]\n")
            for element in settledata:
                if int(element[0]) in np.array(qmatomlist).astype(int):
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


def make_new_top(
    top, molfindlist, mollist, mollength, qmatomlist, includelist, outname
):
    new_mollist = mollist
    molcount = 0
    curr_offset = 0
    extra_count = 0
    for molecule in mollist:
        curr_topname = molfindlist[molcount][1]
        if int(molecule[1]) != 1:
            non_qm_mols = 0
            for j in range(0, int(molecule[1])):
                for atom in qmatomlist:
                    if int(curr_offset) + int(mollength[molcount]) >= int(atom) and int(
                        atom
                    ) > int(curr_offset):
                        # The molecule contains QM atoms and is not unique
                        new_topname = molfindlist[molcount][1]
                        newnew_topname = new_topname + str(extra_count)
                        found = True
                        while found == True:
                            found = False
                            for i in range(0, len(molfindlist)):
                                if molfindlist[i][1] == newnew_topname:
                                    extra_count += 1
                                    newnew_topname = new_topname + str(extra_count)
                                    found = True
                                    break
                        new_mollist.append(
                            write_new_itp(
                                curr_topname,
                                newnew_topname,
                                qmatomlist,
                                curr_offset,
                                mollength[molcount],
                                includelist,
                                mollist,
                            )
                        )
                    else:
                        non_qm_mols += 1
                curr_offset += int(mollength[molcount])
        else:
            # the molecule is unique
            for atom in qmatomlist:
                if int(curr_offset) + int(mollength[molcount]) >= int(atom) and int(
                    atom
                ) > int(curr_offset):
                    # the unique molecule contains a qm atom
                    print("Ding")
                else:
                    print("Dong")
        for i in range(0, int(molecule[1])):
            curr_offset += mollength[molcount]
        molcount += 1


def make_exclude_index(outname, qmatomlist, inner, outer):
    if outer != []:
        ref = qmatomlist + inner + outer # Added by Nicola
        ref = list(set(ref)) # Added by Nicola
        ref.sort() # Added by Nicola
        memory_dict = reindexing.reindexing_memory(ref, outer)
        for i in range(len(qmatomlist)): # Added by Nicola
            qmatomlist[i] = memory_dict[qmatomlist[i]]
        for i in range(len(inner)): # Added by Nicola
            inner[i] = memory_dict[inner[i]]
    with open(outname, "w") as ofile:
        count = 0
        ofile.write("[ QM ]\n")
        count = 0
        for element in qmatomlist:
            count += 1
            ofile.write(str(int(element)) + " ")
            if count == 15:
                ofile.write("\n")
                count = 0
        ofile.write("\n")
        if not inner == []:
            ofile.write("[ INNER ]\n")
            count=0
            for element in inner:
                count += 1
                ofile.write(str(int(element)) + " ")
                if count == 15:
                    ofile.write("\n")
                    count = 0
            ofile.write("\n")
        #if not outer == []:
        #    ofile.write("[ OUTER ]\n")
        #    count=0
        #    for element in outer:
        #            count+=1
        #            ofile.write(str(int(element)) + " ")
        #            if count==15:
        #                    ofile.write("\n")
        #                    count=0
        #    ofile.write("\n")


def generate_top_listsonly(
    top,
    qmatomlist,
    outname,
    flaglist,
    logfile,
    pathinfo,
    inner=[],
    outer=[],
):
    
    #qmatomlist = prep_pcf.read_qmatom_list(qmatoms)
    #inneratomlist = prep_pcf.read_inner_list(inner)
    #outeratomlist = prep_pcf.read_outer_list(outer)
    mollist = make_pcf.readmols(top)
    includelist = make_pcf.getincludelist(top, pathinfo)
    molfindlist = get_molfindlist(mollist, top, includelist)
    make_large_top(
        top,
        molfindlist,
        mollist,
        qmatomlist,
        includelist,
        outname,
        flaglist,
        logfile,
        inner,
        outer
    )

'''
def generate_top(sysargs):
    basedir = os.path.dirname(os.path.abspath(__file__))
    top = sysargs[1]
    qmatoms = sysargs[2]
    outname = sysargs[3]
    flaglist = sysargs[4]
    qmatomlist = prep_pcf.read_qmatom_list(qmatoms)
    mollist = make_pcf.readmols(top)
    includelist = make_pcf.getincludelist(top, pathinfo)
    molfindlist = get_molfindlist(mollist, top, includelist)
    make_large_top(
        top,
        molfindlist,
        mollist,
        qmatomlist,
        includelist,
        outname,
        flaglist,
        logfile,
    )


if __name__ == "__main__":
    import sys

    generate_top(sys.argv)
'''