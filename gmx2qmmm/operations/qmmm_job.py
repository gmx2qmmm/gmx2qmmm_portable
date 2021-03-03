# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will read prepared fields and files for a qmmm job based on gmx.
# will not work alone!

__author__ = "jangoetze"
__date__ = "$06-Feb-2018 12:45:17$"

import math
import os
import re
import subprocess
import copy

import numpy as np
import sqlite3

from gmx2qmmm._helper import logger, _flatten, stepper
from gmx2qmmm._helper import get_linkatoms_ang, make_xyzq
from gmx2qmmm.operations import expansion_check as rot
from gmx2qmmm.operations import nma_stuff
from gmx2qmmm.operations import nma_3N_6dof as nma
from gmx2qmmm.operations import hes_xyz_g09RevD_01_fchk as write_hess
from gmx2qmmm.operations import bmat
from gmx2qmmm.pointcharges import generate_pcf_from_top as make_pcf
from gmx2qmmm.pointcharges import prepare_pcf_for_shift as prep_pcf
from gmx2qmmm.pointcharges import generate_charge_shift as final_pcf

#Ultilites
def get_m2charges(xyzq, m1list, m2list):
    m2charges = []
    count = 0
    for element in m1list:
        m2chargeline = []
        for i in range(0, len(m2list[count])):
            m2chargeline.append(float(xyzq[int(m2list[count][i]) - 1][3]))
        count += 1
        m2charges.append(m2chargeline)
    return m2charges

def read_pcffile(pcffile):
    pcf = []
    with open(pcffile) as ifile:
        for line in ifile:
            match = re.search(r"^QM", line, flags=re.MULTILINE)
            if match:
                pcf.append(["QM"])
                continue
            match = re.search(
                r"^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
            )
            if match:
                pcf.append(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                        float(match.group(4)),
                    ]
                )
                continue
            match = re.search(r"^$end", line, flags=re.MULTILINE)
            if match:
                break
    return pcf

def read_pcf_self(qmfile):
    pcf_self = 0.0
    with open(qmfile + ".log") as ifile:
        for line in ifile:
            match = re.search(
                r"^\s+Self\s+energy\s+of\s+the\s+charges\s+=\s+([-]*\d+\.\d+)\s+a\.u\.",
                line,
                flags=re.MULTILINE,
            )
            if match:
                pcf_self = float(match.group(1))
                break
    return pcf_self

def remove_inactive(total_force, active):
    new_total_force = []
    for i in range(0, len(total_force)):
        if (i + 1) in np.array(active).astype(int):
            new_total_force.append(total_force[i])
        else:
            new_total_force.append([0.0, 0.0, 0.0])
    return new_total_force

def make_clean_force(total_force):

    clean_force = []
    for element in np.array(total_force):
        forceline = []
        for entry in np.array(element):
            forceline.append(float(entry))
        clean_force.append(forceline)
    return clean_force

def get_full_coords_angstrom(gro):
    fullcoords = []
    with open(gro) as ifile:
        count = 0
        for line in ifile:
            count += 1
            if count == 4:
                break
        count = 0
        for line in ifile:
            count += 1
            match = re.search(
                r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if not match:
                break
            full_coord_line = [
                float(match.group(5)) * 10.0,
                float(match.group(6)) * 10.0,
                float(match.group(7)) * 10.0,
            ]
            fullcoords.append(full_coord_line)
    return fullcoords

def get_atoms(qmmmtop, logfile):
    atoms = []
    mass_map = {
        "H": "1.008",
        "He": "4.0026",
        "Li": "6.94",
        "Be": "9.0122",
        "B": "10.81",
        "C": "12.011",
        "N": "14.007",
        "O": "15.999",
        "F": "18.998",
        "Ne": "20.180",
        "Na": "22.990",
        "Mg": "24.305",
        "Al": "26.982",
        "Si": "28.085",
        "P": "30.974",
        "S": "32.06",
        "Cl": "35.45",
        "Ar": "39.948",
        "K": "39.098",
        "Ca": "40.0784",
        "Sc": "44.956",
        "Ti": "47.867",
        "V": "50.942",
        "Cr": "51.996",
        "Mn": "54.938",
        "Fe": "55.8452",
        "Co": "58.933",
        "Ni": "58.693",
        "Cu": "63.5463",
        "Zn": "65.382",
        "Ga": "69.723",
        "Ge": "72.6308",
        "As": "74.922",
        "Se": "78.9718",
        "Br": "79.904",
        "Kr": "83.7982",
        "Rb": "85.468",
        "Sr": "87.62",
        "Y": "88.906",
        "Zr": "91.2242",
        "Nb": "92.906",
        "Mo": "95.95",
        "Tc": "98.906254721",
        "Ru": "101.072",
        "Rh": "102.91",
        "Pd": "106.42",
        "Ag": "107.87",
        "Cd": "112.41",
        "In": "114.82",
        "Sn": "118.71",
        "Sb": "121.76",
        "Te": "127.603",
        "I": "126.90",
        "Xe": "131.29",
        "Cs": "132.91",
        "Ba": "137.33",
        "La": "138.91",
        "Ce": "140.12",
        "Pr": "140.91",
        "Nd": "144.24",
        "Pm": "144.9127493",
        "Sm": "150.362",
        "Eu": "151.96",
        "Gd": "157.253",
        "Tb": "158.93",
        "Dy": "162.50",
        "Ho": "164.93",
        "Er": "167.26",
        "Tm": "168.93",
        "Yb": "173.05",
        "Lu": "174.97",
        "Hf": "178.492",
        "Ta": "180.95",
        "W": "183.84",
        "Re": "186.21",
        "Os": "190.233",
        "Ir": "192.22",
        "Pt": "195.08",
        "Au": "196.97",
        "Hg": "200.59",
        "Tl": "204.38",
        "Pb": "207.2",
        "Bi": "208.98",
        "Po": "208.982430420",
        "At": "209.9871488",
        "Rn": "222.017577725",
        "Fr": "223.019735926",
        "Ra": "226.025409825",
        "Ac": "227.027752126",
        "Th": "232.04",
        "Pa": "231.04",
        "U": "238.03",
        "Np": "237.04817342",
        "Pu": "244.0642045",
        "Am": "243.061381125",
        "Cm": "247.0703545",
        "Bk": "247.0703076",
        "Cf": "251.0795875",
        "Es": "252.082985",
        "Fm": "257.0951067",
        "Md": "258.0984315",
        "No": "259.1010311",
        "Lr": "266.1198356",
        "Rf": "267.1217962",
        "Db": "268.1256757",
        "Sg": "269.1286339",
        "Bh": "270.1333631",
        "Hs": "277.1519058",
        "Mt": "278.1563168",
        "Ds": "281.1645159",
        "Rg": "282.1691272",
        "Cn": "285.177126",
        "Nh": "286.1822172",
        "Fl": "289.190426",
        "Mc": "289.1936389",
        "Lv": "293.204496",
        "Ts": "294.2104674",
        "Og": "295.2162469",
    }
    name_map = {value: key for key, value in mass_map.items()}

    with open(qmmmtop) as ifile:
        for line in ifile:
            match = re.search(r"\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
            if match:
                break
        for line in ifile:
            match = re.search(r"\[\s+atoms\s*\]", line, flags=re.MULTILINE)
            if match:
                break
        for line in ifile:
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
                    if massdiff > 0.01:
                        logger(
                            logfile,
                            str(
                                "Found a mass of "
                                + str(atommass)
                                + " for atom type "
                                + str(atomtype)
                                + ' (identified as atom name "'
                                + str(foundname)
                                + '"), which is more than 0.01 different from the expected mass of '
                                + str(testmass)
                                + ". Atom index was "
                                + str(match.group(1))
                                + ". This has no effect except unless the atom was identified wrongly or dynamics are intended. Clean your ffnonbonded.itp to avoid these messages!\n"
                            ),
                        )
                    atoms.append(foundname)
                else:
                    logger(
                        logfile,
                        str(
                            "Atom type "
                            + str(atomtype)
                            + " could not be translated to a regular atom name. Exiting. Last line:\n"
                        ),
                    )
                    logger(logfile, line)
                    exit(1)
    return atoms

def get_nbradius(gro):

    fullcoords = np.array(get_full_coords_nm(gro))
    mindist = fullcoords[0]
    maxdist = fullcoords[1]
    for element in fullcoords:
        for i in range(0, 3):
            if float(mindist[i]) > float(element[i]):
                mindist[i] = element[i]
            if float(maxdist[i]) < float(element[i]):
                maxdist[i] = element[i]
    maxcoords = np.array(maxdist) - np.array(mindist)
    return np.linalg.norm(maxcoords)

def get_full_coords_nm(gro):  # read g96
    fullcoords = []
    with open(gro) as ifile:
        count = 0
        for line in ifile:
            if count == 4:
                break
            count += 1
        for line in ifile:
            match = re.search(
                r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if not match:
                break
            full_coord_line = [
                float(match.group(5)),
                float(match.group(6)),
                float(match.group(7)),
            ]
            fullcoords.append(full_coord_line)
    return fullcoords

def get_approx_hessian(xyz, old_xyz, grad, old_grad, old_hess, logfile):
    s = xyz - old_xyz
    if s.shape[1] != 1:
        s = s.reshape(3 * len(s), 1)  # reshape s so it can be used in dot products
    g = grad - old_grad
    if g.shape[1] != 1:
        g = g.reshape(3 * len(g), 1)  # reshape g so it can be used in dot products
    # Broyden-Fletcher-Goldfarb-Shanno
    # Define a few values and matrices first for convenience
    mat = old_hess.dot(s)
    factor1 = g.T.dot(s)
    factor2 = s.T.dot(mat)
    idmat = np.eye(len(g))

    # update formula
    new_hess = old_hess + g.dot(g.T) / factor1 - mat.dot(mat.T) / factor2
    # moreover, we calculate eigenvalues as they are indicative of the curvature of the current PES
    eigvals, eigvecs = np.linalg.eig(new_hess)

    # Check BFGS condition
    if factor1 > 0:
        logger(logfile, "BFGS condition fulfilled.\n")
        WARN = False
    else:
        logger(logfile, "BFGS condition not fulfilled! We keep the old Hessian.\n")
        WARN = True
        new_hess = old_hess

    return new_hess, eigvals.min(), WARN

#qm & mm program ultis
def make_g16_inp(qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    gro = qmmmInputs.gro
    qmmmtop = qmmmInputs.top
    #qmmmtop = qmmmInputs.qmmmtop
    qmatomlist = qmmmInputs.qmatomlist
    pcffile = qmmmInputs.pcffile
    curr_step = qmmmInputs.qmmmparams.curr_step
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    nmaflag = qmmmInputs.nmaflag

    method = qmmmInputs.qmparams.method
    basis = qmmmInputs.qmparams.basis
    charge = qmmmInputs.qmparams.charge
    multi = qmmmInputs.qmparams.multi
    cores = qmmmInputs.qmparams.cores
    memory = qmmmInputs.qmparams.memory
    extra = qmmmInputs.qmparams.extra
    
    insert = ""
    oldinsert = ""

    gaufile = str(jobname + insert + ".gjf")
    chkfile = str(jobname + insert + ".chk")
    oldchkfile = str(jobname + oldinsert + ".chk")

    if nmaflag == 1:
        oldchkfile = str(jobname + ".chk")

    with open(gaufile, "w") as ofile:
        fullcoords = get_full_coords_angstrom(gro)
        atoms = get_atoms(qmmmtop, logfile)
        ofile.write("%NPROCSHARED=" + str(cores) + "\n")
        ofile.write("%MEM=" + str(memory) + "MB\n")
        ofile.write("%CHK=" + chkfile + "\n")
        if int(curr_step) != 0 or nmaflag == 1:
            ofile.write("%OLDCHK=" + oldchkfile + "\n")
        ofile.write("#P " + str(method))
        if str(basis) != "NONE":
            ofile.write("/" + str(basis))
        if str(extra) != "NONE":
            ofile.write(" " + str(extra))
        if int(curr_step) != 0 or nmaflag == 1:
            ofile.write(" guess=read")
        ofile.write(
            " nosymm gfinput gfprint force charge guess=huckel punch=derivatives iop(3/33=1) prop(field,read) pop=esp\n"
        )
        ofile.write(
            "\nQMMM Calc QM part\n\n"
            + str(int(charge))
            + " "
            + str(int(multi))
            + "\n"
        )
        count = 0

        for element in fullcoords:
            if int(count + 1) in np.array(qmatomlist).astype(int):
                ofile.write(
                    "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                        str(atoms[count]),
                        float(element[0]),
                        float(element[1]),
                        float(element[2]),
                    )
                )
            count += 1
        for element in linkatoms:
            # links are already in angstrom
            ofile.write(
                "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                    str("H"), float(element[0]), float(element[1]), float(element[2])
                )
            )
        ofile.write("\n")

        with open(pcffile) as ifile:
            for line in ifile:
                match = re.search(
                    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                )
                if match:
                    ofile.write(
                        "{:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                            float(match.group(1)),
                            float(match.group(2)),
                            float(match.group(3)),
                            float(match.group(4)),
                        )
                    )
        ofile.write("\n")

        with open(pcffile) as ifile:
            for line in ifile:
                match = re.search(
                    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                )
                if match:
                    ofile.write(
                        "{:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                            float(match.group(1)),
                            float(match.group(2)),
                            float(match.group(3)),
                        )
                    )
        ofile.write("\n\n")

    return gaufile

def make_gmx_inp(qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    gro = qmmmInputs.gro
    qmmmtop = qmmmInputs.qmmmtop
    qmatomlist = qmmmInputs.qmatomlist
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    basedir = qmmmInputs.basedir
    prefix =  qmmmInputs.pathparams.gmxpath + qmmmInputs.pathparams.gmxcmd

    insert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step)))
    mdpname = str(jobname + ".mdp")
    groname = str(jobname + ".boxlarge.g96")
    ndxname = str(qmmmtop + ".ndx")
    tprname = str(jobname + insert + ".tpr")
    nbradius = get_nbradius(gro)
    write_mdp(mdpname, nbradius)
    update_gro_box(gro, groname, nbradius, logfile)
    subprocess.call(
        [
            prefix,
            "grompp",
            "-p",
            str(qmmmtop),
            "-c",
            str(groname),
            "-n",
            str(ndxname),
            "-f",
            str(mdpname),
            "-o",
            str(tprname),
            "-backup",
            "no",
        ]
    )
    subprocess.call(["rm", "mdout.mdp"])
    return tprname

def write_mdp(mdpname, nbradius):
    with open(mdpname, "w") as ofile:
        ofile.write(
            "title               =  Yo\ncpp                 =  /usr/bin/cpp\nconstraints         =  none\nintegrator          =  md\ndt                  =  0.001 ; ps !\nnsteps              =  1\nnstcomm             =  0\nnstxout             =  1\nnstvout             =  1\nnstfout             =   1\nnstlog              =  1\nnstenergy           =  1\nnstlist             =  1\nns_type             =  grid\nrlist               =  "
        )
        ofile.write(str(float(nbradius)))
        ofile.write(
            "\ncutoff-scheme = group\ncoulombtype    =  cut-off\nrcoulomb            =  "
        )
        ofile.write(str(float(nbradius)))
        ofile.write("\nrvdw                =  ")
        ofile.write(str(float(nbradius)))
        ofile.write(
            "\nTcoupl              =  no\nenergygrps          =  QM\nenergygrp-excl = QM QM\nPcoupl              =  no\ngen_vel             =  no\n"
        )

def update_gro_box(gro, groname, nbradius, logfile):
    with open(groname, "w") as ofile:
        with open(gro) as ifile:
            logger(
                logfile,
                str(
                    "Finding a larger .gro box size to avoid problems with .mdp input..."
                ),
            )
            for line in ifile:
                ofile.write(line)
                match = re.search(r"^BOX\s*\n", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\d*\.\d+)\s+(\d*\.\d+)\s+(\d*\.\d+)",
                            line,
                            flags=re.MULTILINE,
                        )
                        if not match:
                            logger(
                                logfile,
                                "\n\nError: In "
                                + str(gro)
                                + " box vectors were expected but not found. Exiting. Line was:\n",
                            )
                            logger(logfile, line)
                            exit(1)
                        else:
                            bv = [
                                float(match.group(1)) + 10.0 * nbradius,
                                float(match.group(2)) + 10.0 * nbradius,
                                float(match.group(3)) + 10.0 * nbradius,
                            ]
                            ofile.write(
                                " {:>15.9f} {:>15.9f} {:>15.9f}\nEND\n".format(
                                    float(bv[0]), float(bv[1]), float(bv[2])
                                )
                            )
                        break
                    break
    logger(logfile, str("Done.\n"))

#Propagated
#initstep = stepsize
def propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, curr_step, logfile, scan_flag=False):
    dispvec = []
    maxforce = 0.0
    clean_force = make_clean_force(total_force)
    old_clean_force = make_clean_force(last_forces)
    maxatom = -1
    maxcoord = -1
    check_force = list(_flatten(clean_force))
    #search index of max info
    for i in range(0, int(len(check_force) / 3)):
        for j in range(0, 3):
            if abs(float(check_force[i * 3 + j])) > abs(maxforce):
                maxforce = float(check_force[i * 3 + j])
                maxatom = i
                maxcoord = j
    logger(
        logfile,
        str(
            "Maximum force is "
            + str(float(maxforce))
            + " a.u. at coord "
            + str(int(maxatom) + 1)
            + "/"
            + str(int(maxcoord) + 1)
            + ".\n"
        ),
    )
    #All use steep for first propagation
    if curr_step <= 1:
        logger(logfile, "Propagate with steepest descent for first step...\n")
        for element in clean_force:
            dispvec.append(
                [
                    float(element[0]) * float(stepsize) / abs(float(maxforce)),
                    float(element[1]) * float(stepsize) / abs(float(maxforce)),
                    float(element[2]) * float(stepsize) / abs(float(maxforce)),
                ]
            )
    else:
        if propagator == "STEEP":
            logger(logfile, "Propagate with steepest descent...\n")
            for element in clean_force:
                dispvec.append(
                    [
                        float(element[0]) * float(stepsize) / abs(float(maxforce)),
                        float(element[1]) * float(stepsize) / abs(float(maxforce)),
                        float(element[2]) * float(stepsize) / abs(float(maxforce)),
                    ]
                )
            corr_length = np.array(total_force)

        elif propagator == "CONJGRAD" :
            logger(logfile, "Propagate with conjugate gradient...\n")
            # Fletcher-Reeves
            _flattened = list(_flatten(clean_force))
            corr_fac = np.array(_flattened).dot(np.array(_flattened))
            _flattened = list(_flatten(old_clean_force))
            corr_fac /= np.array(_flattened).dot(np.array(_flattened))

            #2222
            corr_length = np.array(total_force)
            if curr_step > 1:   #0 step for the inital SP 
                corr_length = np.array(total_force) + corr_fac * corr_length

            for i in range(len(total_force)):
                dispvec.append(
                    [
                        float(stepsize) * corr_length[i][0],
                        float(stepsize) * corr_length[i][1],
                        float(stepsize) * corr_length[i][2],
                    ]
                )

            logger(
                logfile,
                str(
                    "Effective step at maximum force coord is "
                    + str(float(dispvec[maxatom][maxcoord]))
                    + " a.u.\n"
                ),
            )
        elif propagator == "BFGS":
            logger(logfile, "Propagate with BFGS method...\n")
            coords = np.array(new_xyzq)[:, 0:3]
            old_coords = np.array(xyzq)[:, 0:3]
            old_hessian = np.loadtxt("bfgs_hessian.txt")
            gradient = np.array(total_force)
            old_gradient = np.array(last_forces)
            hessian, hesseig, warning_flag = get_approx_hessian(
                coords, old_coords, gradient, old_gradient, old_hessian, logfile
            )
            np.savetxt("bfgs_hessian.txt", hessian)

            coords = coords.reshape(
                3 * len(coords), 1
            )  # reshape coords to use in dot products
            gradient = gradient.reshape(
                3 * len(gradient), 1
            )  # reshape grad to use in dot products

            direc = -np.linalg.inv(hessian).dot(gradient)
            direc = direc.reshape(int(len(coords) / 3), 3)

            for i in range(len(total_force)):
                dispvec.append(
                    [
                        float(stepsize) * direc[i][0],
                        float(stepsize) * direc[i][1],
                        float(stepsize) * direc[i][2],
                    ]
                )
    

    return dispvec

def make_g96_inp(dispvec, gro, new_gro, logfile):

    with open(new_gro, "w") as ofile:
        with open(gro) as ifile:
            counter = 0
            for line in ifile:
                ofile.write(line)
                counter += 1
                if counter == 4:
                    break
            counter = 0
            for line in ifile:
                match = re.search(
                    r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    ofile.write(line)
                    logger(
                        logfile,
                        str(
                            "Successfully wrote "
                            + str(int(counter))
                            + " atoms to new g96 file.\n"
                        ),
                    )
                    break
                else:
                    dispx = dispvec[counter][0] * 0.052917721
                    dispy = dispvec[counter][1] * 0.052917721
                    dispz = dispvec[counter][2] * 0.052917721
                    ofile.write(
                        str(match.group(1))
                        + " "
                        + str(match.group(2))
                        + " "
                        + str(match.group(3))
                        + " "
                        + str(match.group(4))
                        + " {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                            float(match.group(5)) + float(dispx),
                            float(match.group(6)) + float(dispy),
                            float(match.group(7)) + float(dispz),
                        )
                    )
                    counter += 1
            for line in ifile:
                ofile.write(line) 
    
    ifile.close()
    ofile.close()
    logger(logfile, str("Made new coordinates in file " + str(new_gro) + ".\n"))   

def qmmm_prep(qmmmInputs):
    gro = qmmmInputs.gro
    top = qmmmInputs.top
    jobname = qmmmInputs.qmmmparams.jobname
    qmatomlist = qmmmInputs.qmatomlist
    connlist = qmmmInputs.connlist
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    basedir = qmmmInputs.basedir
    gmxtop = qmmmInputs.pathparams.gmxtop 
    charge = qmmmInputs.qmparams.charge
    geo = make_pcf.readg96(gro)
    curr_step = qmmmInputs.qmmmparams.curr_step

    logger(logfile, "List of molecules...\n")
    # 6277
    mollist = make_pcf.readmols(top)
    logger(logfile, "Done.\n")
    logger(logfile, "Reading charges...\n")
    chargevec = []
    for element in mollist:
        chargevec.extend(make_pcf.readcharges(element, top, gmxtop))
    logger(logfile, "Done.\n")
    new_xyzq = make_xyzq(geo, chargevec)
    logger(logfile, str("Made new xyzq matrix.\n"))
    logger(
        logfile,
        "Preparing the point charge field for a numerically optimized charge shift...\n",
    )
    (
        qmcoordlist,
        m1list,
        m2list,
        updated_chargelist,
    ) = prep_pcf.prepare_pcf_for_shift_fieldsonly(
        new_xyzq, qmatomlist, charge, connlist
    )
    logger(logfile, "Done.\n")
    new_links = get_linkatoms_ang(
        new_xyzq, qmatomlist, m1list, connlist, linkatoms
    )
    logger(logfile, str("Updated positions of link atoms.\n"))
    filename = jobname
    if curr_step > 0:
        filename += "." + str(int(curr_step))
    if not os.path.isfile(str(filename + ".pointcharges")):
        logger(logfile, "Shifting...\n")
        final_pcf.generate_charge_shift_fieldsonly(
            updated_chargelist, m1list, qmcoordlist, m2list, filename, basedir
        )
        logger(logfile, str("Made new PCF file.\n"))
    else:
        logger(
            logfile,
            "NOTE: Shifting omitted due to "
            + str(filename + ".pointcharges")
            + " being an existing file!\n",
        )
    logger(logfile, "Done.\n")

    #Update new output
    qmmmInputs.xyzq = new_xyzq
    qmmmInputs.m1list = m1list
    qmmmInputs.m2list = m2list
    qmmmInputs.linkatoms = new_links


    return qmmmInputs

#Database
def databasecorrection(energy_or_force, cut, dist, qmmmInputs):
    
    forcefield = qmmmInputs.mmparams.fField
    method = qmmmInputs.qmparams.method
    basisset = qmmmInputs.qmparams.basis
    fit = qmmmInputs.qmmmparams.databasefit
    basedir = qmmmInputs.basedir
    logfile = qmmmInputs.logfile

    conn = sqlite3.connect(basedir + "/correction_database/database.sqlite")

    # check if method exist in database
    method_set = conn.cursor()
    method_set.execute(
        "SELECT * FROM "
        + cut
        + ' WHERE forcefield="'
        + forcefield
        + '" AND method="'
        + method
        + '" AND basisset ="'
        + basisset
        + '"'
    )
    method_set_value = method_set.fetchall()
    if len(method_set_value) == 0:
        cut = "aminoacid_CACB"
        forcefield = "amberGS"
        method = "CAM-B3LYP"
        basisset = "6-31++G**"
        logger(
            logfile,
            "Unexisted method in correction database, changing to default correction method...\n",
        )

    c = conn.cursor()
    c.execute(
        "SELECT * FROM "
        + cut
        + ' WHERE forcefield="'
        + forcefield
        + '" AND method="'
        + method
        + '" AND basisset ="'
        + basisset
        + '"'
    )
    db_values = c.fetchall()[0]

    conn.close()
    returnvalue = 0
    if len(db_values) > 0:
        if energy_or_force == "ENERGY":
            if fit == "POLY":
                returnvalue = (
                    db_values[5] * dist * dist * dist
                    + db_values[6] * dist * dist
                    + db_values[7] * dist
                    + db_values[8]
                )
            elif fit == "MORSE":
                returnvalue = (
                    db_values[9]
                    * (
                        np.exp(-2 * db_values[10] * (dist - db_values[11]))
                        - 2 * np.exp(-db_values[10] * (dist - db_values[11]))
                    )
                    + db_values[12]
                )
            elif fit == "NO":
                returnvalue = 0
                logger(logfile, "No energy correction.\n")
        elif energy_or_force == "FORCES":
            if fit == "POLY":
                returnvalue = (
                    db_values[13] * dist * dist * dist
                    + db_values[14] * dist * dist
                    + db_values[15] * dist
                    + db_values[16]
                )
            elif fit == "MORSE":
                returnvalue = (
                    db_values[17]
                    * (
                        np.exp(-2 * db_values[18] * (dist - db_values[19]))
                        - 2 * np.exp(-db_values[18] * (dist - db_values[19]))
                    )
                    + db_values[20]
                )
            elif fit == "NO":
                returnvalue = 0
                logger(logfile, "No force correction.\n")
    return returnvalue

#Run program command
def run_g16(qmfile, qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    g16cmd = qmmmInputs.pathparams.g16cmd
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    #2222
    insert = ""
    #if int(curr_step) > 0:
    #    insert = str("." + str(int(curr_step)))


    if not os.path.isfile(str(qmfile) + ".log"):
        logger(logfile, "Running G16 file.\n")
        subprocess.call([g16cmd, str(qmfile)])
        logname = qmfile[:-3]
        logname += "log"
        subprocess.call(["mv", logname, str(jobname + insert + ".gjf.log")])
        subprocess.call(["mv", "fort.7", str(jobname + insert + ".fort.7")])
        logger(logfile, "G16 Done.\n")
    else:
        logger(
            logfile,
            "NOTE: Using existing G16 files, skipping calculation for this step.\n",
        )
    if not os.path.isfile(jobname + insert + ".fort.7"):
        if not os.path.isfile("fort.7"):
            logger(
                logfile,
                "No fort.7 file was created by the last Gaussian run! Exiting.\n",
            )
            exit(1)
        subprocess.call(["mv", "fort.7", str(jobname + insert + ".fort.7")])
        logger(
            logfile,
            "WARNING: Had to rename fort.7 file but not the log file. MAKE SURE THAT THE FORT.7 FILE FITS TO THE LOG FILE!\n",
        )

def run_gmx(mmfile, qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    prefix =  qmmmInputs.pathparams.gmxpath + qmmmInputs.pathparams.gmxcmd
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile

    insert = ""
    if int(curr_step) != 0:
        insert = str("." + str(int(curr_step)))

    logger(logfile, "Running Gromacs file.\n")
    trrname = str(jobname + insert + ".trr")
    xtcname = str(jobname + insert + ".xtc")
    outname = str(jobname + insert + ".out.gro")
    gmxlogname = str(jobname + insert + ".gmx.log")
    edrname = str(jobname + insert + ".edr")
    subprocess.call(
        [
            prefix,
            "mdrun",
            "-s",
            mmfile,
            "-o",
            trrname,
            "-c",
            outname,
            "-x",
            xtcname,
            "-g",
            gmxlogname,
            "-e",
            edrname,
            "-backup",
            "no",
        ]
    )
    subprocess.call(["rm", outname])

    return edrname

# Write output file
    #old output
def old_write_output(energies, total_force, curr_step):
    qmenergy, mmenergy, linkcorrenergy, total_energy = energies
    energy_file = "oenergy.txt"
    
    if curr_step == 0:
        file_flag = "w"
    else:
        
        file_flag = "a+"
    
    oenergy = open(energy_file, file_flag)
    oenergy.write("Step: %d\n" % curr_step)
    oenergy.write("QM energy: %.8f\n" % qmenergy)
    oenergy.write("MM energy: %.8f\n" % mmenergy)
    oenergy.write("Link energy: %.8f\n" % linkcorrenergy)
    oenergy.write("Total energy: %.8f\n"% total_energy)
    oenergy.write("--------------------------------------------\n")
    oenergy.close()

    oforce = open("oforce.txt", file_flag)
    for i in range(len(total_force)):
        oforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
        )
    oforce.write('\n')
    oforce.close()
    #test output

def write_test(qmmmInputs, curr_step):
    qmenergy, mmenergy, linkcorrenergy, total_energy = qmmmInputs.energies
    total_force = qmmmInputs.forces
    #qmenergy, mmenergy, linkcorrenergy, total_energy = energies
    energy_file = "testE.txt"
    
    if curr_step == 0:
        file_flag = "w"
    else:
        
        file_flag = "a+"
    
    oenergy = open(energy_file, file_flag)
    oenergy.write("Step: %d\n" % curr_step)
    oenergy.write("QM energy: %.8f\n" % qmenergy)
    oenergy.write("MM energy: %.8f\n" % mmenergy)
    oenergy.write("Link energy: %.8f\n" % linkcorrenergy)
    oenergy.write("Total energy: %.8f\n"% total_energy)
    oenergy.write("--------------------------------------------\n")
    oenergy.close()

    oforce = open("testF.txt", file_flag)
    for i in range(len(total_force)):
        oforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
        )
    oforce.write('\n')
    oforce.close()

def write_output(energies, total_force, curr_step, energy_file, forces_file):
    qmenergy, mmenergy, linkcorrenergy, total_energy = energies
    #energy_file = "oenergy.txt"
    if curr_step == 0:
        file_flag = "w"
    else:
        file_flag = "a+"

    oenergy = open(energy_file, file_flag)
    
    if curr_step == 0:
        oenergy.write("Step\tQM\t\tMM\t\tLink\t\tTotal\n")
        oenergy.write("%d\t%f\t%f\t%f\t%f\n" % (curr_step, qmenergy, mmenergy, linkcorrenergy, total_energy))
    else:
        
        oenergy.write("%d\t%f\t%f\t%f\t%f\n" % (curr_step, qmenergy, mmenergy, linkcorrenergy, total_energy))
    oenergy.close()
    #forces_file = "oforces.txt"
    oforce = open(forces_file, file_flag)
    oforce.write("Step%d\n"%curr_step)
    for i in range(len(total_force)):
        oforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
        )
    oforce.write("\n")
    oforce.close()

# Read output file
def read_oenergy(step):
    oenergy = np.loadtxt('oenergy.txt',dtype=float,skiprows=1,usecols=(1,2,3,4))
    return oenergy[step]

def read_oforces(step):
    file_info = open('oforces.txt', 'r')
    step_index = []
    for i, line in enumerate(file_info.readlines()): 
        if 'Step' in line:
            step_index.append(i)
    steps = len(step_index)
    with open('oforces.txt','r') as file_object:
        line = file_object.readlines()
    for i in range(len(step_index)):
        if i < len(step_index)-1:
            outlines = step_index[i+1] - step_index[i]
        else:
            outlines = len(line) - step_index[i] - 1
    forces = []
    for i in range(outlines):
        force = line[step_index[step]+i+1].split()[1:]
        for j in range(3):
            force[j] = float(force[j])
        forces.append(force)
    return np.array(forces)

# Archive & remove
def archive(jobname, curr_step):
    if curr_step == 0 :
        insert = '' 
    else:
        insert = str("." + str(curr_step))

    output_filename = str(jobname + insert + ".tar")
    archive_files = str(jobname + insert + "*")  
    subprocess.call("tar -czf %s %s"%(output_filename, archive_files), shell=True)

def remove_files(jobname, curr_step, tar=False):
    if curr_step == 0 :
        insert = '' 
    else:
        insert = str("." + str(curr_step))

    extend = np.array([".g96",".boxlarge.g96",".chk",".edr",".edr.xvg",".fort.7",".gjf",".gjf.log",
        ".gmx.log", ".mdp", ".pointcharges", ".qmmm.top",".qmmm.top.ndx", ".tpr", ".trr", ".xvg"])

    for i in range(len(extend)):
        subprocess.call("rm -f %s"%str(jobname + insert + extend[i]), shell=True)

    if tar:
        subprocess.call("rm -f %s"%str(jobname + insert + '.tar'), shell=True)

#Energy
def get_qmenergy(qmfile, qmmmInputs):
    qmprog = qmmmInputs.qmparams.program
    extra_string = qmmmInputs.qmparams.extra
    pcffile = qmmmInputs.pcffile
    logfile = qmmmInputs.logfile
    basedir = qmmmInputs.basedir

    logger(logfile, "Extracting QM energy.\n")
    qmenergy = 0.0
    qm_corrdata = []
    if str(qmprog) == "G16":
        with open(str(qmfile + ".log")) as ifile:
            for line in ifile:
                match = []
                match2 = []
                match2 = re.search(
                    r"\sTD[=(\s]", extra_string.upper(), flags=re.MULTILINE
                )
                if not match2:
                    match2 = re.search(
                        r"^TD[=(\s]", extra_string.upper(), flags=re.MULTILINE
                    )
                if not match2:
                    match2 = re.search(
                        r"\sTD$", extra_string.upper(), flags=re.MULTILINE
                    )
                if not match2:
                    match2 = re.search(
                        r"^TD$", extra_string.upper(), flags=re.MULTILINE
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
                    logger(logfile, "Obtaining charge self-interaction...\n")
                    pcf_self_pot = read_pcf_self(qmfile)
                    logger(
                        logfile, "Done: {:>20.10f} a.u.\n".format(float(pcf_self_pot))
                    )
                    # G16 energy needs to be corrected for self potential of PCF
                    qmenergy = float(match.group(1)) - float(pcf_self_pot)
                match = re.search(r"^\s*ESP\s*charges:", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        break
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\d+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:
                            qm_corrdata.append(
                                [
                                    int(match.group(1)),
                                    match.group(2),
                                    float(match.group(3)),
                                ]
                            )
                        else:
                            break
                    break
    logger(logfile, "QM energy is " + str(float(qmenergy)) + " a.u..\n")
    return qmenergy, qm_corrdata

def get_mmenergy(edrname, qmmmInputs):
    prefix =  qmmmInputs.pathparams.gmxpath + qmmmInputs.pathparams.gmxcmd
    logfile = qmmmInputs.logfile

    mmenergy = 0.0
    logger(logfile, "Extracting MM energy.\n")
    p = subprocess.Popen(
        [
            prefix,
            "energy",
            "-f",
            edrname,
            "-o",
            str(edrname + ".xvg"),
            "-backup",
            "no",
        ],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    p.communicate(input=b"11\n\n")
    with open(str(edrname + ".xvg")) as ifile:
        for line in ifile:
            match = re.search(
                r"^    0.000000\s*([-]*\d+.\d+)\n", line, flags=re.MULTILINE
            )
            if match:
                mmenergy = float(match.group(1)) * 0.00038087988
                break
    logger(logfile, "MM energy is " + str(float(mmenergy)) + " a.u..\n")
    return mmenergy

def get_linkenergy_au(qm_corrdata, qmmmInputs):
    xyzq = qmmmInputs.xyzq
    linkcorrlist = qmmmInputs.linkcorrlist
    m1list = qmmmInputs.m1list
    m2list = qmmmInputs.m2list
    q1list = qmmmInputs.q1list
    qmmmtop = qmmmInputs.qmmmtop
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    pcffile = qmmmInputs.pcffile
    qmatomlist = qmmmInputs.qmatomlist

    linkenergy = 0.0
    m2charges = get_m2charges(xyzq, m1list, m2list)
    for element in linkcorrlist:
        z1 = 0.0
        v1 = []
        v2 = []
        if int(element[0]) in np.array(list(_flatten(m2list))).astype(int):
            for i in range(0, len(m2list)):
                for j in range(0, len(m2list[i])):
                    if int(m2list[i][j]) == int(element[0]):
                        z1 = float(m2charges[i][j])
                        v1 = [
                            xyzq[int(element[0]) - 1][0] / 0.52917721,
                            xyzq[int(element[0]) - 1][1] / 0.52917721,
                            xyzq[int(element[0]) - 1][2] / 0.52917721,
                        ]
                        break
                if z1 != 0.0:
                    break
        elif int(element[0]) in np.array(qmatomlist).astype(int):
            for i in range(0, len(qmatomlist)):
                if int(qmatomlist[i]) == int(element[0]):
                    z1 = float(qm_corrdata[i][2])
                    v1 = [
                        xyzq[int(element[0]) - 1][0] / 0.52917721,
                        xyzq[int(element[0]) - 1][1] / 0.52917721,
                        xyzq[int(element[0]) - 1][2] / 0.52917721,
                    ]
                    break
        elif int(element[0]) in np.array(m1list).astype(int):
            for i in range(0, len(m1list)):
                if int(m1list[i]) == int(element[0]):
                    z1 = float(qm_corrdata[i + len(qmatomlist)][2])
                    v1 = [
                        linkatoms[i][0] / 0.52917721,
                        linkatoms[i][1] / 0.52917721,
                        linkatoms[i][2] / 0.52917721,
                    ]
                    break
        else:
            z1 = float(xyzq[int(element[0]) - 1][3])
            v1 = [
                xyzq[int(element[0]) - 1][0] / 0.52917721,
                xyzq[int(element[0]) - 1][1] / 0.52917721,
                xyzq[int(element[0]) - 1][2] / 0.52917721,
            ]
        z2 = 0.0
        if int(element[1]) in _flatten(m2list):
            for i in range(0, len(m2list)):
                for j in range(0, len(m2list[i])):
                    if int(m2list[i][j]) == int(element[1]):
                        z2 = float(m2charges[i][j])
                        v2 = [
                            xyzq[int(element[1]) - 1][0] / 0.52917721,
                            xyzq[int(element[1]) - 1][1] / 0.52917721,
                            xyzq[int(element[1]) - 1][2] / 0.52917721,
                        ]
                        break
                if z2 != 0.0:
                    break
        elif int(element[1]) in np.array(qmatomlist).astype(int):
            for i in range(0, len(qmatomlist)):
                if int(qmatomlist[i]) == int(element[1]):
                    z2 = float(qm_corrdata[i][2])
                    v2 = [
                        xyzq[int(element[1]) - 1][0] / 0.52917721,
                        xyzq[int(element[1]) - 1][1] / 0.52917721,
                        xyzq[int(element[1]) - 1][2] / 0.52917721,
                    ]
                    break
        elif int(element[1]) in np.array(m1list).astype(int):
            for i in range(0, len(m1list)):
                if int(m1list[i]) == int(element[1]):
                    z2 = float(qm_corrdata[i + len(qmatomlist)][2])
                    v2 = [
                        linkatoms[i][0] / 0.52917721,
                        linkatoms[i][1] / 0.52917721,
                        linkatoms[i][2] / 0.52917721,
                    ]
                    break
        else:
            z2 = float(xyzq[int(element[1]) - 1][3])
            v2 = [
                xyzq[int(element[1]) - 1][0] / 0.52917721,
                xyzq[int(element[1]) - 1][1] / 0.52917721,
                xyzq[int(element[1]) - 1][2] / 0.52917721,
            ]
        v12 = np.array(v1) - np.array(v2)
        dist = np.linalg.norm(v12)
        linkenergy += z1 * z2 / dist
    # now also all atoms in the corrdata list with the mod and linkcorr point charges
    # mod first. mod is charge in pcffile minus m2charge
    pcf = read_pcffile(pcffile)
    for i in range(0, len(m2list)):
        for j in range(0, len(m2list[i])):
            curr_mod = []
            for k in range(0, 3):
                curr_mod.append(float(pcf[int(m2list[i][j]) - 1][k]) / 0.52917721)
            curr_mod_charge = (
                float(float(pcf[int(m2list[i][j]) - 1][3])) - m2charges[i][j]
            )
            for k in range(0, len(qmatomlist)):
                v1 = [
                    xyzq[int(qmatomlist[k]) - 1][0] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][1] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                linkenergy += z1 * curr_mod_charge / dist
            for k in range(0, len(linkatoms)):
                v1 = [
                    linkatoms[k][0] / 0.52917721,
                    linkatoms[k][1] / 0.52917721,
                    linkatoms[k][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k + len(qmatomlist)][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                linkenergy += z1 * curr_mod_charge / dist
    # now linkcorr. linkcorr are last m2*2 entries in pcf
    m2count = 0
    linkstart = len(pcf) - 2 * len(list(_flatten(m2list)))
    for i in range(0, len(m2list)):
        for j in range(0, len(m2list[i])):
            curr_mod = []
            for k in range(0, 3):
                curr_mod.append(float(pcf[int(linkstart) + m2count][k]) / 0.52917721)
            curr_mod_charge = float(float(pcf[int(linkstart) + m2count][3]))
            m2count += 1
            for k in range(0, len(qmatomlist)):
                v1 = [
                    xyzq[int(qmatomlist[k]) - 1][0] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][1] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                linkenergy += z1 * curr_mod_charge / dist
            for k in range(0, len(linkatoms)):
                v1 = [
                    linkatoms[k][0] / 0.52917721,
                    linkatoms[k][1] / 0.52917721,
                    linkatoms[k][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k + len(qmatomlist)][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                linkenergy += z1 * curr_mod_charge / dist
    # now, add the correction of energy for the link atoms. currently only C-C bond cuts supported.
    for i in range(0, len(linkatoms)):
        v1 = [
            linkatoms[i][0] / 0.52917721,
            linkatoms[i][1] / 0.52917721,
            linkatoms[i][2] / 0.52917721,
        ]
        _flattened = list(_flatten(q1list))
        v2 = [
            xyzq[int(_flattened[i]) - 1][0] / 0.52917721,
            xyzq[int(_flattened[i]) - 1][1] / 0.52917721,
            xyzq[int(_flattened[i]) - 1][2] / 0.52917721,
        ]
        v12 = np.array(v2) - np.array(v1)
        dist = np.linalg.norm(v12)
    dist = dist * 0.7409471631
    energycorr = databasecorrection("ENERGY", "aminoacid_CACB", dist, qmmmInputs)
    linkenergy -= energycorr
    # sign inverted due to correction convention (subtracting)
    return linkenergy

def get_energy(qmfile, edrname, qmmmInputs):
    logfile = qmmmInputs.logfile
    qmenergy, qm_corrdata = get_qmenergy(qmfile, qmmmInputs)
    mmenergy = get_mmenergy(str(edrname), qmmmInputs)
    linkcorrenergy = get_linkenergy_au(qm_corrdata, qmmmInputs)
    basis = qmmmInputs.qmparams.basis
    #qmenergy -= linkcorrenergy
    methodstring = str(qmmmInputs.qmparams.method)
    if basis != "NONE":
        methodstring += str("/" + str(basis))
    logger(
        logfile,
        str(
            "Single point energy done. QM/MM energy is {:>20.10f} (QM, link atom corrected ".format(
                float(qmenergy-linkcorrenergy)
            )
            + methodstring
            + ") + {:>20.10f} (MM) = {:>20.10f} (a.u.)\n".format(
                float(mmenergy), float(qmenergy-linkcorrenergy) + float(mmenergy)
            )
        ),
    )

    return qmenergy, mmenergy, linkcorrenergy, qm_corrdata

#Forces
def get_qmforces_au(qmmmInputs):
    qmatomlist = qmmmInputs.qmatomlist
    m1list = qmmmInputs.m1list
    qmmmtop = qmmmInputs.qmmmtop
    qmprogram = qmmmInputs.qmparams.program
    jobname = qmmmInputs.qmmmparams.jobname
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile

    qmforces = []
    qmonlyforcelist = []
    pcf_grad = []

    if qmprogram == "G16":
        insert = ""
        #2222
        if (int(curr_step) > 1) and (qmmmInputs.qmmmparams.jobname == 'OPT'):
            insert = str("." + str(curr_step-1))
        
        qmlogfile = str(jobname + insert + ".gjf.log")
        fortfile = str(jobname + insert + ".fort.7")

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
    with open(qmmmtop) as ifile:
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
                                    not in np.array(qmatomlist).astype(int)
                                )
                                and (int(match.group(1)) not in np.array(m1list).astype(int))
                            ):
                                curr_charge = float(match.group(7))
                                qmforces.append(
                                    [
                                        pcf_grad[count][0] * curr_charge,
                                        pcf_grad[count][1] * curr_charge,
                                        pcf_grad[count][2] * curr_charge,
                                    ]
                                )
                                count += 1
                            elif match and int(match.group(1)) in np.array(
                                qmatomlist
                            ).astype(int):
                                qmforces.append(qmonlyforcelist[qmcount])
                                qmcount += 1
                            elif match and int(match.group(1)) in np.array(m1list).astype(
                                int
                            ):
                                qmforces.append(
                                    qmonlyforcelist[m1count + len(qmatomlist)]
                                )
                                m1count += 1
                            match = re.search(r"^\s*\n", line, flags=re.MULTILINE)
                            if match:
                                break
                        break
                break
    return qmforces

def get_mmforces_au(qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    curr_step = qmmmInputs.qmmmparams.curr_step
    prefix =  qmmmInputs.pathparams.gmxpath + qmmmInputs.pathparams.gmxcmd
    logfile = qmmmInputs.logfile

    mmforces = []
    insert = ""
    if int(curr_step) != 0:
        insert = str("." + str(curr_step))
    trrname = str(jobname + insert + ".trr")
    tprname = str(jobname + insert + ".tpr")
    xvgname = str(jobname + insert + ".xvg")
    p = subprocess.Popen(
        [
            prefix,
            "traj",
            "-fp",
            "-f",
            trrname,
            "-s",
            tprname,
            "-of",
            xvgname,
            "-xvg",
            "none",
            "-backup",
            "no",
        ],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    p.communicate(input=b"0\n")
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
                    mmforces.append(mmforceline)
                    mmforceline = []
            break  # read only one line

    return mmforces

def get_linkforces_au(qm_corrdata,qmmmInputs):
    linkcorrlist = qmmmInputs.linkcorrlist
    qmatomlist = qmmmInputs.qmatomlist
    xyzq = qmmmInputs.xyzq
    pcffile = qmmmInputs.pcffile
    m1list = qmmmInputs.m1list
    m2list = qmmmInputs.m2list
    q1list = qmmmInputs.q1list
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile

    linkforces = []
    # Force Coulomb: z1*z2*(distance along coord)/(distance between charges)**3
    for element in xyzq:  # this is just to count an entry for each atom!
        linkforces.append([0.0, 0.0, 0.0])
    m2charges = get_m2charges(xyzq, m1list, m2list)
    for element in linkcorrlist:
        z1 = 0.0
        v1 = []
        v2 = []
        if int(element[0]) in _flatten(m2list):
            for i in range(0, len(m2list)):
                for j in range(0, len(m2list[i])):
                    if int(m2list[i][j]) == int(element[0]):
                        z1 = float(m2charges[i][j])
                        v1 = [
                            xyzq[int(element[0]) - 1][0] / 0.52917721,
                            xyzq[int(element[0]) - 1][1] / 0.52917721,
                            xyzq[int(element[0]) - 1][2] / 0.52917721,
                        ]
                        break
                if z1 != 0.0:
                    break
        elif int(element[0]) in np.array(qmatomlist).astype(int):
            for i in range(0, len(qmatomlist)):
                if int(qmatomlist[i]) == int(element[0]):
                    z1 = float(qm_corrdata[i][2])
                    v1 = [
                        xyzq[int(element[0]) - 1][0] / 0.52917721,
                        xyzq[int(element[0]) - 1][1] / 0.52917721,
                        xyzq[int(element[0]) - 1][2] / 0.52917721,
                    ]
                    break
        elif int(element[0]) in np.array(m1list).astype(int):
            for i in range(0, len(m1list)):
                if int(m1list[i]) == int(element[0]):
                    z1 = float(qm_corrdata[i + len(qmatomlist)][2])
                    v1 = [
                        linkatoms[i][0] / 0.52917721,
                        linkatoms[i][1] / 0.52917721,
                        linkatoms[i][2] / 0.52917721,
                    ]
                    break
        else:
            z1 = float(xyzq[int(element[0]) - 1][3])
            v1 = [
                xyzq[int(element[0]) - 1][0] / 0.52917721,
                xyzq[int(element[0]) - 1][1] / 0.52917721,
                xyzq[int(element[0]) - 1][2] / 0.52917721,
            ]
        z2 = 0.0
        if int(element[1]) in _flatten(m2list):
            for i in range(0, len(m2list)):
                for j in range(0, len(m2list[i])):
                    if int(m2list[i][j]) == int(element[1]):
                        z2 = float(m2charges[i][j])
                        v2 = [
                            xyzq[int(element[1]) - 1][0] / 0.52917721,
                            xyzq[int(element[1]) - 1][1] / 0.52917721,
                            xyzq[int(element[1]) - 1][2] / 0.52917721,
                        ]
                        break
                if z2 != 0.0:
                    break
        elif int(element[1]) in np.array(qmatomlist).astype(int):
            for i in range(0, len(qmatomlist)):
                if int(qmatomlist[i]) == int(element[1]):
                    z2 = float(qm_corrdata[i][2])
                    v2 = [
                        xyzq[int(element[1]) - 1][0] / 0.52917721,
                        xyzq[int(element[1]) - 1][1] / 0.52917721,
                        xyzq[int(element[1]) - 1][2] / 0.52917721,
                    ]
                    break
        elif int(element[1]) in np.array(m1list).astype(int):
            for i in range(0, len(m1list)):
                if int(m1list[i]) == int(element[1]):
                    z2 = float(qm_corrdata[i + len(qmatomlist)][2])
                    v2 = [
                        linkatoms[i][0] / 0.52917721,
                        linkatoms[i][1] / 0.52917721,
                        linkatoms[i][2] / 0.52917721,
                    ]
                    break
        else:
            z2 = float(xyzq[int(element[1]) - 1][3])
            v2 = [
                xyzq[int(element[1]) - 1][0] / 0.52917721,
                xyzq[int(element[1]) - 1][1] / 0.52917721,
                xyzq[int(element[1]) - 1][2] / 0.52917721,
            ]
        v12 = np.array(v1) - np.array(v2)
        dist = np.linalg.norm(v12)

        for i in range(0, 3):
            linkforces[int(element[0]) - 1][i] += (
                z1 * z2 * v12[i] / (dist * dist * dist)
            )
            linkforces[int(element[1]) - 1][i] -= (
                z1 * z2 * v12[i] / (dist * dist * dist)
            )
    # now also all atoms in the corrdata list with the mod and linkcorr point charges
    # mod first. mod is charge in pcffile minus m2charge
    pcf = read_pcffile(pcffile)
    for i in range(0, len(m2list)):
        for j in range(0, len(m2list[i])):
            curr_mod = []
            for k in range(0, 3):
                curr_mod.append(float(pcf[int(m2list[i][j]) - 1][k]) / 0.52917721)
            curr_mod_charge = (
                float(float(pcf[int(m2list[i][j]) - 1][3])) - m2charges[i][j]
            )
            for k in range(0, len(qmatomlist)):
                v1 = [
                    xyzq[int(qmatomlist[k]) - 1][0] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][1] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                for l in range(0, 3):
                    linkforces[int(qmatomlist[k]) - 1][l] += (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
                    linkforces[int(m2list[i][j]) - 1][l] -= (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
            for k in range(0, len(linkatoms)):
                v1 = [
                    linkatoms[k][0] / 0.52917721,
                    linkatoms[k][1] / 0.52917721,
                    linkatoms[k][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k + len(qmatomlist)][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                for l in range(0, 3):
                    linkforces[int(m1list[k]) - 1][l] += (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
                    linkforces[int(m2list[i][j]) - 1][l] -= (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
    m2count = 0
    linkstart = len(pcf) - 2 * len(list(_flatten(m2list)))
    for i in range(0, len(m2list)):
        for j in range(0, len(m2list[i]) * 2):
            curr_mod = []
            for k in range(0, 3):
                curr_mod.append(float(pcf[int(linkstart) + m2count][k]) / 0.52917721)
            curr_mod_charge = float(float(pcf[int(linkstart) + m2count][3]))
            m2count += 1
            for k in range(0, len(qmatomlist)):
                v1 = [
                    xyzq[int(qmatomlist[k]) - 1][0] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][1] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                for l in range(0, 3):
                    linkforces[int(qmatomlist[k]) - 1][l] += (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
            for k in range(0, len(linkatoms)):
                v1 = [
                    linkatoms[k][0] / 0.52917721,
                    linkatoms[k][1] / 0.52917721,
                    linkatoms[k][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k + len(qmatomlist)][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                for l in range(0, 3):
                    linkforces[int(m1list[k]) - 1][l] += (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
    for i in range(0, len(linkatoms)):
        v1 = [
            linkatoms[i][0] / 0.52917721,
            linkatoms[i][1] / 0.52917721,
            linkatoms[i][2] / 0.52917721,
        ]
        _flattened = list(_flatten(q1list))
        v2 = [
            xyzq[int(_flattened[i]) - 1][0] / 0.52917721,
            xyzq[int(_flattened[i]) - 1][1] / 0.52917721,
            xyzq[int(_flattened[i]) - 1][2] / 0.52917721,
        ]
        v12 = np.array(v2) - np.array(v1)
        dist = np.linalg.norm(v12) / 0.71290813568205
        u_v12 = rot.uvec(v12)
        dist = dist * 0.5282272551
        forcecorr = databasecorrection("FORCES", "aminoacid_CACB", dist, qmmmInputs)
        for j in range(0, 3):
            linkforces[int(_flattened[i]) - 1][j] += -u_v12[j] * forcecorr * 0.5
            linkforces[int(m1list[i]) - 1][j] += u_v12[j] * forcecorr * 0.5
    return linkforces

def read_forces(qm_corrdata,qmmmInputs):
    logfile = qmmmInputs.logfile
    active = qmmmInputs.active
    qmforces = []
    mmforces = []
    qmforces = get_qmforces_au(qmmmInputs)
    logger(logfile, str("QM forces read.\n"))
    mmforces = get_mmforces_au(qmmmInputs)
    logger(logfile, str("MM forces read.\n"))
    linkcorrforces = get_linkforces_au(qm_corrdata, qmmmInputs)
    logger(logfile, str("Forces for link atom correction read.\n"))
    total_force = np.array(qmforces) + np.array(mmforces) - np.array(linkcorrforces)
    logger(logfile, str("Total forces obtained.\n"))
    if (qmmmInputs.qmmmparams.jobtype != "SINGLEPOINT") and (len(active) != 0):
        total_force=remove_inactive(total_force,active) #in SP case I don't have inactive atoms
        logger(logfile, str("Deleted forces of inactive atoms.\n"))
    return total_force

#update input
def make_opt_step(qmmmInputs):
    return 0

#scan
def scan_dispvec(xyzq, stepsize, scan_atoms):
    dispvec = np.zeros((len(xyzq),3))
    coords = np.array(xyzq)[:, 0:3]
    #Bond length
    if len(scan_atoms) == 2 :
        a,b = scan_atoms
        coord_A, coord_B = (coords[int(a-1)], coords[int(b-1)])# array index -> -1
        line_vec = bmat.length(coord_A, coord_B)*stepsize
        dispvec[int(b-1)] += line_vec
    
    #Angle
    elif len(scan_atoms) == 3 :
        angle_vec = bmat.angle(coords[scan_atoms[0]], coords[scan_atoms[1]], coords[scan_atoms[2]])

    #Dihedral
    elif len(scan_atoms) == 4 :
        dihedral = bmat.angle(coords[scan_atoms[0]], coords[scan_atoms[1]], coords[scan_atoms[2]], coords[scan_atoms[3]])

    return dispvec

#Job
def perform_sp(qmmmInputs):
    jobtype = qmmmInputs.qmmmparams.jobtype
    qmprogram = qmmmInputs.qmparams.program
    logfile = qmmmInputs.logfile
    curr_step = qmmmInputs.qmmmparams.curr_step

    logger(logfile, "Computing a single point.\n")
    logger(logfile, "Preparing QM and MM inputs:\n ")
    if jobtype == "OPT":
        logger(logfile, "-----OPT %d---------------.\n"%curr_step)

    if qmprogram == "G16":
        #qm
        qmfile = make_g16_inp(qmmmInputs)
        logger(logfile, "Gaussian input file, %s, is ready.\n"%qmfile)
        print('-----Run g16file:%s---------------\n'%qmfile)
        run_g16(qmfile, qmmmInputs)

        #mm
        mmfile = make_gmx_inp(qmmmInputs)
        logger(logfile, "Gromacs input file, %s, is ready.\n"%mmfile)      
        edrname = run_gmx(mmfile, qmmmInputs)
        
        logger(logfile, str("Reading energies.\n"))
        qmenergy, mmenergy, linkcorrenergy, qm_corrdata = get_energy(qmfile, edrname, qmmmInputs)
        total_energy = qmenergy + mmenergy - linkcorrenergy
        energies = (qmenergy, mmenergy, linkcorrenergy, total_energy)
        logger(logfile, str("Reading forces.\n"))
        total_force = read_forces(qm_corrdata,qmmmInputs)

        qmmmInputs.energies = energies
        qmmmInputs.forces = total_force

    elif qmprogram == "TM":
        logger(logfile,"Turbomole is not avalible currently.\n")
        exit(0)

    elif qmprogram == "ORCA":
        logger(logfile,"ORCA is not avalible currently.\n")
        exit(0)
    else:
        logger(logfile,"Unidentified QM program. Please check your -qm intput.\n")
        exit(0)

    # write a total force file
    #if jobtype == "SINGLEPOINT" or curr_step == 0:
    if jobtype == "SINGLEPOINT":
        logger(logfile, "Writing SP output.\n")
        write_output(energies, total_force, curr_step, "oenergy.txt", "oforces.txt")

    elif jobtype == "OPT":
        logger(logfile, "Writing OPT output.\n")

    elif jobtype == "NMA":
        logger(logfile, "Writing NMA output.\n")

    elif jobtype == "SCAN":
        logger(logfile, "Writing SCAN output.\n")
 
def perform_opt(qmmmInputs):

    #define done
    STEPLIMIT = 0
    FTHRESH = 1
    STEPSIZE = 2 

    logfile = qmmmInputs.logfile
    propagator = qmmmInputs.qmmmparams.propagater
    curr_step = qmmmInputs.qmmmparams.curr_step
    maxcycle = qmmmInputs.qmmmparams.maxcycle
    f_thresh = qmmmInputs.qmmmparams.f_thresh
    optlastonly = qmmmInputs.qmmmparams.optlastonly
    stepsize = qmmmInputs.qmmmparams.initstep
    jobname = qmmmInputs.qmmmparams.jobname
    jobtype = qmmmInputs.qmmmparams.jobtype
    qmprog = qmmmInputs.qmparams.program

    gro = qmmmInputs.gro
    xyzq = qmmmInputs.xyzq
    new_xyzq = []
    count = 0 #for each calculation counting for archiving

    if curr_step == 0:
        #First calculation 
        perform_sp(qmmmInputs)
        old_qmmmInputs = copy.deepcopy(qmmmInputs)
        
        if jobtype == "SCAN" :
            write_output(qmmmInputs.energies, qmmmInputs.forces, curr_step, "oenergy_%s.txt"%jobname, "oforces_%s.txt"%jobname)
        else:
            write_output(qmmmInputs.energies, qmmmInputs.forces, curr_step, "oenergy.txt", "oforces.txt")
        
        #Store 1st result
        curr_energy = qmmmInputs.energies[-1] #total energy
        total_force = qmmmInputs.forces

        # init BFGS hessian: initial guess for identity hessian
        if propagator == "BFGS":
            init_hessian = np.eye(3 * len(xyzq))
            np.savetxt("bfgs_hessian.txt", init_hessian)
    else:
        curr_energy = read_oenergy(curr_step) #total energy
        total_force = read_oforces(curr_step)
    
    last_forces = []
    clean_force = make_clean_force(total_force)
    maxforce = 0.0
    done = STEPLIMIT
    improved = True

    for element in _flatten(clean_force):
        if abs(float(element)) > abs(maxforce):
            maxforce = float(element)
    if abs(maxforce) < float(f_thresh):
        logger(logfile,"Max force (%f) below threshold (%f) Finishing.\n"%(maxforce,f_thresh))
        done = FTHRESH    
    else:
        logger(logfile, "Max force not below threshold. Continuing.\n")

        ############## start optimization loop ##############
        while not done and count <= maxcycle:
            #Store previous result
            old_qmmmInputs = copy.deepcopy(qmmmInputs)

            
            #Add step first (if not improve, reject -1)
            qmmmInputs.qmmmparams.curr_step += 1
            curr_step = qmmmInputs.qmmmparams.curr_step
            
            #Prepare new input            
            new_gro = str(jobname + "." + str(curr_step) + ".g96")
            new_pcffile = str(jobname + "." + str(curr_step) + ".pointcharges")
            qmmmInputs.gro = new_gro
            qmmmInputs.pcffile = new_pcffile
            
            if jobtype == "SCAN" :
                dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, curr_step, logfile, True)
            else :
                dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, curr_step, logfile)
            
            make_g96_inp(dispvec, gro, new_gro, logfile)
            qmmm_prep(qmmmInputs)

            #Run SP
            print("00000000000",qmmmInputs.pcffile)
            perform_sp(qmmmInputs) #make g16 & gmx input

            #Store new result
            curr_energy = qmmmInputs.energies[-1] #total energy
            last_energy = old_qmmmInputs.energies[-1]
            total_force = qmmmInputs.forces
            last_forces = old_qmmmInputs.forces

            ############## start energy check & update stepsize ############## 
            if curr_energy > last_energy:
                logger(logfile, "Rejected one optimization step due to energy increasing. Trying again, with smaller step.\n")
                improved = False
                stepsize *= 0.2

                #remove files                
                insert = str("." + str(curr_step))
                trrname = str(jobname + insert + ".trr")
                tprname = str(jobname + insert + ".tpr")
                gmxlogname = str(jobname + insert + ".gmx.log")
                edrname = str(jobname + insert + ".edr")
                xvgname = str(jobname + insert + ".edr.xvg")
                g16name = str(jobname + insert + ".gjf.log")
                fortname = str(jobname + insert + ".fort.7")
                subprocess.call("rm %s %s %s %s %s %s %s"%(trrname,tprname,gmxlogname,edrname,xvgname,g16name,fortname), shell=True)
                

            else:
                stepsize *= 1.2
                improved = True

                #archive remove previous
                logger(logfile, "Archive previous files\n")
                archive(jobname, curr_step-1) #archive previous


            ############## end energy check & update stepsize ##############
            
            ############## start force check ##############
            clean_force = make_clean_force(total_force)
            maxforce = 0.0
            
            for element in _flatten(clean_force):
                if abs(float(element)) > abs(maxforce):
                    maxforce = float(element)
            if abs(maxforce) < float(f_thresh):
                logger(logfile,"Max force (%f) below threshold (%f) Finishing.\n"%(maxforce,f_thresh))
                done = FTHRESH
                break
            ############## end force check ##############

            ############## stepsize check ##############

            if float(stepsize) < 1e-6: #0.000001 a.u.
                done = STEPSIZE
                logger(
                    logfile,
                        ("Step became lower than 0.000001 a.u., optimization is considered done for now. " + 
                         "This is the best we can do unless reaching unacceptable numerical noise levels.\n"),
                )
                break
            ############## end stepsize check ##############

            if improved:
                #Update for BFGS              
                xyzq = old_qmmmInputs.xyzq
                new_xyzq = qmmmInputs.xyzq
                #Store previous 
                old_qmmmInputs = copy.deepcopy(qmmmInputs)
                count += 1

                if jobtype == "SCAN" :
                    write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step, "oenergy_%s.txt"%jobname, "oforces_%s.txt"%jobname)
                else:
                    write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step, "oenergy.txt", "oforces.txt")
            
            else:
                #rejected and use prvious
                qmmmInputs.qmmmparams.curr_step -= 1
                qmmmInputs = copy.deepcopy(old_qmmmInputs)


        ############## end optimization loop ##############
    #Remain first/last result 
    if qmmmInputs.qmmmparams.optlastonly == "YES":
        logger(logfile, "Remain last result. Remove other results.\n")
        #remove first
        remove_files(jobname, 0, True)
        for steps in range(count-1):
            filename = jobname+('.%d*'%(steps+1))
            subprocess.call("rm %s"%(filename), shell=True)
    else:
        for steps in range(count):
            remove_files(jobname,steps)

    # Remove the last but failed with criteria files         
    subprocess.call("rm %s"%str(jobname + '.' + str(count+1) + '*'), shell=True)
    
    # opt status 
    if done == STEPLIMIT:
        logger(logfile, "Optimization canceled due to step limit.\n")
    elif done == FTHRESH:
        logger(logfile, "Optimization finished due to force threshold.\n")
    elif done == STEPSIZE:
        logger(logfile, "Optimization finished due to step size.\n")
   
def perform_scan(qmmmInputs):
    logfile = qmmmInputs.logfile
    qmmmInputs.qmmmparams.propagator = "BFGS"
    r_array, a_array, d_array = qmmmInputs.qmmmparams.scan
    
    #Check scan input
    if (len(r_array)+len(a_array)+len(d_array)) == 0:
        logger(logfile, "No scan coordinates are found. Please check the scan file input.\n")
    else:
        logger(logfile, 'Read %d scan coordinates.\n'%(len(r_array)+len(a_array)+len(d_array)))

    #Single scan
    #Bond scan
    if len(r_array) > 0:
        r_shape = r_array.shape
        print("r_shape", r_shape, '\n')
        subprocess.call("mkdir scanR", shell=True)
        origin_qmmmInputs = copy.deepcopy(qmmmInputs)
        # Bond length scan loop
        for i in range(r_shape[0]):
            scan_atoms = (r_array[i][0], r_array[i][1])
            stepsize = r_array[0][-2]
            steps = r_array[0][-1]
            logger(logfile, "Scanning bond length between atom%d-%d\n"%scan_atoms)
            logger(logfile, "Scanned stepsize: %.3f, scan %d steps\n"%(stepsize, steps))
            
            direc = "scanR/R%d-%d"%scan_atoms
            subprocess.call("mkdir scanR/R%d-%d"%scan_atoms, shell=True)
            origin_gro = qmmmInputs.gro
            # Bond length step loop
            steps=2
            for j in range(steps):
                qmmmInputs = copy.deepcopy(origin_qmmmInputs)                
                logger(logfile, "Start scanned step%d...\n"%(j+1))

                logger(logfile, "Create new scanned geometry with increament %.3f\n"%(stepsize*(j+1)) )
                dispvec = scan_dispvec(qmmmInputs.xyzq, (stepsize*(j+1)), scan_atoms)
                new_gro = "scanR%d-%d"%scan_atoms + "_%d.g96"%(j+1)
                make_g96_inp(dispvec, origin_gro, new_gro, logfile)
                
                qmmmInputs.gro = new_gro
                qmmmInputs.qmmmparams.jobname = "scanR%d-%d"%scan_atoms + "_%d"%(j+1)
                qmmmInputs.qmmmparams.curr_step = 0
                
                logger(logfile, ("Run optimization of scanR%d-%d"%scan_atoms + "_%d\n"%(j+1)) )
                perform_opt(qmmmInputs)

                logger(logfile, "End scanned step%d...\n"%(j+1))
           

    #Angle scan
    if len(a_array) > 0:
        a_shape = a_array.shape
        print("a_shape", a_shape, '\n')
        subprocess.call("mkdir scanA", shell=True)
        for i in range(a_shape[0]):
            logger(logfile, "Scanning angle between atom%d-%d-%d\n"%(a_array[i][0], a_array[i][1], a_array[i][2]))
            logger(logfile, "Scanned stepsize: %f, scan %d steps\n"%(a_array[0][-2], a_array[0][-1]))
            subprocess.call("mkdir scanA/A%d-%d-%d"%(a_array[i][0], a_array[i][1], a_array[i][2]), shell=True)

    #Dihedral angle scan
    if len(d_array) > 0:
        d_shape = d_array.shape
        print("d_shape", d_shape, '\n')
        subprocess.call("mkdir scanD", shell=True)
        for i in range(d_shape[0]):
            logger(logfile, "Scanning dihedral angle between atom%d-%d-%d-%d\n"%(d_array[i][0], d_array[i][1], d_array[i][2], d_array[i][3]))
            logger(logfile, "Scanned stepsize: %f, scan %d steps\n"%(d_array[0][-2], d_array[0][-1]))
            subprocess.call("mkdir scanD/D%d-%d-%d-%d"%(d_array[i][0], d_array[i][1], d_array[i][2], d_array[i][3]), shell=True)

def perform_nma(qmmmInputs):
    logfile = qmmmInputs.qmmmparams.logfile
    basedir = qmmmparams.basedir
    active = qmmmInputs.active
    jobname = qmmmInputs.qmmmparams.jobname

    logger(logfile, "------This will be a numerical) normal mode analysis.------\n")
    logger(
        logfile,
        "Generating a numerical Hessian for the active region using a displacement step of %f a.u.\n"%qmmmInputs.qmmmparams.disp,
    )

    logger(logfile, "Will require %d single point calculations!\n"%(len(active)*6+1))

    perform_sp(qmmmInputs)

    start_energy = qmmmInputs.energies[-1]
    start_forces = qmmmInputs.forces
    start_grad = np.array(start_forces) * -1.0

    hessian_xyz_full = []

    for curr_atom in active:
        grad_deriv_vec = nma_stuff.get_xyz_2nd_deriv(qmmmInputs, curr_atom, start_energy, start_forces)
        hessian_xyz_full.extend(grad_deriv_vec)
    prep_hess = nma_stuff.prepare_hess(hessian_xyz_full, active)

    for i in range(0, len(prep_hess[0]) - 6):
        evals.extend([float(1000.0)])
    
    nma_stuff.write_pseudofchk_file(
        jobname, evals, prep_hess, prep_hess, active, qmmmtop, logfile, xyzq
    )  # using prep_hess as pseudo-nm_matrix since we do not know yet
    
    logger(logfile, "Wrote pseudofchk (G03 format).\n")
    
    write_hess.hes_xyz_fchk(
        str(jobname + ".pseudofchk"), str(jobname + ".hess")
    )
    
    logger(logfile, "Wrote orca format .hess file.\n")
    evals, nm_matrix = nma.nma_3Nminus6dof_asfunction(
        str(jobname + ".hess"), basedir
    )
    print(nma_stuff.log_nma(
        qmmminfo, logfile, evals, nm_matrix, active, qmmmtop, xyzq, prep_hess
    ))

    return 0


if __name__ == "__main__":
    print("This file serves as a library for gmx2qmmm-related functions.")
    print("Do not execute directly. Use gmx2qmmm instead.")
