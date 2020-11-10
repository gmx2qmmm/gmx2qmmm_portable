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

import numpy as np
import sqlite3

from gmx2qmmm._helper import logger, _flatten, stepper
from gmx2qmmm._helper import get_linkatoms_ang, make_xyzq
from gmx2qmmm.operations import expansion_check as rot
from gmx2qmmm.operations import scan as scan_func
from gmx2qmmm.operations import nma_stuff
from gmx2qmmm.operations import nma_3N_6dof as nma
from gmx2qmmm.operations import hes_xyz_g09RevD_01_fchk as write_hess
from gmx2qmmm.pointcharges import generate_pcf_from_top as make_pcf
from gmx2qmmm.pointcharges import prepare_pcf_for_shift as prep_pcf
from gmx2qmmm.pointcharges import generate_charge_shift as final_pcf




def perform_sp(qmmmInputs):
    qmmmInputs.qmmmparams.jobname = stepper(qmmmInputs.qmmmparams.jobname, step)
    jobtype = qmmmInputs.qmmmparams.jobtype
    qmprogram = qmmmInputs.qmparams
    logfile = qmmmInputs.logfile

    logger(logfile, "Computing a single point.\n")
    logger(logfile, "Preparing QM and MM inputs: ")

    if qmprogram == "G16":
        
        qmfile = make_g16_inp(qmmmInputs)
        logger(logfile, "Gaussian input ready.\n")
        mmfile, edrname = make_gmx_inp(qmmmInputs)
        logger(logfile, "Gromacs input ready.\n")

        run_g16(qmfile, qmmmInputs)
        edrname = run_gmx(mmfile, qmmmInputs)

        qmenergy, mmenergy, linkcorrenergy, qm_corrdata = get_energy(edrname,qmmmInputs)
        total_force = read_forces(qmmmInputs, qmenergy, mmenergy, qm_corrdata)

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
    if jobtype == "SINGLEPOINT":
        logger(logfile, "Writing SP output.\n")

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
        oenergy.write(
            "Total energy: %.8f\n"
            % (float(qmenergy) + float(mmenergy) - float(linkcorrenergy))
        )
        oenergy.write("--------------------------------------------\n")
        oenergy.close()
        


        oforce = open("oforce.txt", "w")
        for i in range(len(total_force)):
            oforce.write(
                "%d %.8f %.8f %.8f\n"
                % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
            )
        oforce.close()

    elif jobtype == "OPT":
        logger(logfile, "Writing OPT output.\n")

    elif jobtype == "NMA":
        logger(logfile, "Writing NMA output.\n")

    elif jobtype == "SCAN":
        logger(logfile, "Writing SCAN output.\n")
 

def perform_opt(qmmmInputs):

    count = qmmminfo[7]
    qmenergy, mmenergy, qm_corrdata = get_energy(
        gro,
        qmmminfo[0],
        xyzq,
        connlist,
        qmatomlist,
        m1list,
        m2list,
        q1list,
        qmmmtop,
        qminfo,
        mminfo,
        qmmminfo,
        linkcorrlist,
        flaglist,
        pcffile,
        linkatoms,
        int(count),
        logfile,
        basedir,
        pathinfo,
    )
    last_energy = float(qmenergy) + float(mmenergy)
    done = 0
    maxcycle = int(qmmminfo[3])
    initstep = float(qmmminfo[4])
    new_links = linkatoms
    new_pcffile = pcffile
    new_gro = gro
    new_xyzq = xyzq
    new_qm_corrdata = qm_corrdata
    last_forces = []
    while not done and count <= maxcycle:
        logger(
            logfile, str("-----Optimization cycle " + str(int(count) + 1) + "-----\n")
        )
        jobname = qmmminfo[0]
        if count > 0:
            jobname += "." + str(int(count))
        archivename = str(jobname) + ".tar.gz"
        if os.path.isfile(archivename):
            subprocess.call(["tar", "-xf", archivename])
            subprocess.call(["rm", archivename])
        archive = ["tar", "-cf", str(jobname) + ".tar"]
        files = [
            new_pcffile,
            str(jobname) + ".edr",
            str(jobname) + ".edr.xvg",
            str(jobname) + ".trr",
            str(jobname) + ".xvg",
            str(jobname) + ".gmx.log",
            str(jobname) + ".g96",
            str(jobname) + ".gjf",
            str(jobname) + ".mdp",
            str(jobname) + ".boxlarge.g96",
            str(jobname) + ".tpr",
            str(jobname) + ".gjf.log",
            str(jobname) + ".fort.7",
        ]
        (
            done,
            last_energy,
            new_gro,
            new_pcffile,
            new_xyzq,
            new_links,
            initstep,
            new_qm_corrdata,
            last_forces,
            imporved,
        ) = opt_cycle(
            new_gro,
            top,
            new_xyzq,
            connlist,
            qmatomlist,
            new_qm_corrdata,
            m1list,
            m2list,
            q1list,
            qmmmtop,
            qminfo,
            mminfo,
            qmmminfo,
            linkcorrlist,
            flaglist,
            new_pcffile,
            new_links,
            count,
            last_energy,
            last_forces,
            initstep,
            active,
            logfile,
            basedir,
            pathinfo,
        )
        if not imporved:
            logger(logfile, str("No imporvement. Optimization stopped.\n"))
            done = 2
            break
        archive.extend(files)
        subprocess.call(archive)
        subprocess.call(["gzip", str(jobname) + ".tar"])
        delete = ["rm"]
        delete.extend(files)
        subprocess.call(delete)

        energy_file = "optenergy.txt"
        force_file = "optforce.txt"
        rm_file = "orm.txt"
        rmname = qmmminfo[0] + "." + str(int(count))
        if count == 0:
            file_flag = "w"
        elif count > 0:
            file_flag = "a+"
            if qmmminfo[9] == "YES":
                orm = open(rm_file, file_flag)
                orm.write("%s\n" % rmname)
                orm.close()
        if count >= 0:
            oenergy = open(energy_file, file_flag)
            oenergy.write("%d %.8f \n" % (count + 1, last_energy))
            oenergy.close()

            oforce = open(force_file, file_flag)
            oforce.write("Step %d \n" % (count + 1))
            for i in range(len(last_forces)):
                oforce.write(
                    "%d %.8f %.8f %.8f\n"
                    % ((i + 1), last_forces[i][0], last_forces[i][1], last_forces[i][2])
                )
            oforce.close()

        if not done:
            if qmmminfo[9] == "YES":
                subprocess.call(["rm", rmname + ".chk"])
                subprocess.call(["rm", rmname + ".tar.gz"])

        count += 1

    subprocess.call(["rm", qmmminfo[0] + "." + str(int(count + 1)) + ".chk"])
    subprocess.call(["rm", qmmminfo[0] + "." + str(int(count + 1)) + ".gjf"])

    if done == 0:
        logger(logfile, "Optimization canceled due to step limit.\n")
    elif done == 1:
        logger(logfile, "Optimization finished due to energy threshold.\n")
    elif done == 2:
        logger(logfile, "Optimization finished due to step size.\n")
    g96_to_gro(
        str(qmmminfo[0] + "." + str(count) + ".g96"),
        str(qmmminfo[0] + ".opt.gro"),
        logfile,
    )
    logger(logfile, "Final geometry written to " + str(qmmminfo[0]) + ".opt.gro.\n")

    optforce = open("oforce.txt", "w")
    for i in range(len(last_forces)):
        optforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), last_forces[i][0], last_forces[i][1], last_forces[i][2])
        )
    optforce.close()


def opt_cycle(
    gro,
    top,
    xyzq,
    connlist,
    qmatomlist,
    qm_corrdata,
    m1list,
    m2list,
    q1list,
    qmmmtop,
    qminfo,
    mminfo,
    qmmminfo,
    linkcorrlist,
    flaglist,
    pcffile,
    linkatoms,
    curr_step,
    last_energy,
    last_forces,
    initstep,
    active,
    logfile,
    basedir,
    pathinfo,
):

    done = 0
    f_thresh = float(qmmminfo[2])
    jobname = str(qmmminfo[0])
    total_force = read_forces(
        qmatomlist,
        m1list,
        qmmmtop,
        qminfo,
        qmmminfo[0],
        curr_step,
        logfile,
        linkcorrlist,
        xyzq,
        pcffile,
        qm_corrdata,
        m2list,
        q1list,
        linkatoms,
        active,
        basedir,
        mminfo,
        qmmminfo,
        pathinfo,
    )
    clean_force = make_clean_force(total_force)
    maxforce = 0.0
    for element in _flatten(clean_force):
        if abs(float(element)) > abs(maxforce):
            maxforce = float(element)
    if abs(maxforce) < float(f_thresh):
        logger(
            logfile,
            str(
                "Max force "
                + str(maxforce)
                + " below threshold ("
                + str(f_thresh)
                + "). Finishing.\n"
            ),
        )
        done = 1
        return (
            done,
            last_energy,
            gro,
            pcffile,
            xyzq,
            linkatoms,
            initstep,
            qm_corrdata,
            clean_force,
        )
    else:
        logger(logfile, str("Max force not below threshold. Continuing.\n"))

    (
        qmenergy,
        mmenergy,
        new_qm_corrdata,
        new_gro,
        new_pcffile,
        new_xyzq,
        new_links,
        new_initstep,
        imporved,
    ) = make_opt_step(
        jobname,
        xyzq,
        connlist,
        qmmminfo[5],
        gro,
        top,
        total_force,
        initstep,
        qmatomlist,
        qm_corrdata,
        m1list,
        m2list,
        q1list,
        qmmmtop,
        qminfo,
        mminfo,
        qmmminfo,
        linkcorrlist,
        flaglist,
        pcffile,
        linkatoms,
        curr_step,
        last_energy,
        last_forces,
        logfile,
        basedir,
        pathinfo,
    )
    curr_energy = float(qmenergy) + float(mmenergy)
    if float(new_initstep) < 0.000001:
        done = 2
        logger(
            logfile,
            str(
                "Step became lower than 0.000001 a.u., optimization is considered done for now. This is the best we can do unless reaching unacceptable numerical noise levels.\n"
            ),
        )
    if float(curr_energy) > float(last_energy):
        done = 1
        logger(
            logfile,
            str(
                "Energy did not drop! Exiting optimizer, this might indicate an error!\n"
            ),
        )
    return (
        done,
        curr_energy,
        new_gro,
        new_pcffile,
        new_xyzq,
        new_links,
        new_initstep,
        new_qm_corrdata,
        clean_force,
        imporved,
    )

def perform_scan(qmmmInputs):
    # SP
    jobname = stepper(qmmminfo[0], step)
    qmenergy, mmenergy, qm_corrdata = get_energy(
        gro,
        jobname,
        xyzq,
        connlist,
        qmatomlist,
        m1list,
        m2list,
        q1list,
        qmmmtop,
        qminfo,
        mminfo,
        qmmminfo,
        linkcorrlist,
        flaglist,
        pcffile,
        linkatoms,
        int(0),
        logfile,
        basedir,
        pathinfo,
    )
    total_force = read_forces(
        qmatomlist,
        m1list,
        qmmmtop,
        qminfo,
        jobname,
        int(0),
        logfile,
        linkcorrlist,
        xyzq,
        pcffile,
        qm_corrdata,
        m2list,
        q1list,
        linkatoms,
        active,
        basedir,
        mminfo,
        qmmminfo,
        pathinfo,
    )

    scan_data = scan_func.load_scan("scan.txt")

    scan_cycle(
        scan_data,
        gro,
        top,
        xyzq,
        connlist,
        qmatomlist,
        qm_corrdata,
        m1list,
        m2list,
        q1list,
        qmmmtop,
        qminfo,
        mminfo,
        qmmminfo,
        linkcorrlist,
        flaglist,
        pcffile,
        linkatoms,
        curr_step,
        last_energy,
        last_forces,
        initstep,
        active,
        logfile,
        basedir,
        pathinfo,
    )

def scan_cycle(
    scan_data,
    gro,
    top,
    xyzq,
    connlist,
    qmatomlist,
    qm_corrdata,
    m1list,
    m2list,
    q1list,
    qmmmtop,
    qminfo,
    mminfo,
    qmmminfo,
    linkcorrlist,
    flaglist,
    pcffile,
    linkatoms,
    curr_step,
    last_energy,
    last_forces,
    initstep,
    active,
    logfile,
    basedir,
    pathinfo,
):

    # already run sp once

    jobname = str(qmmminfo[0])
    total_force = read_forces(
        qmatomlist,
        m1list,
        qmmmtop,
        qminfo,
        qmmminfo[0],
        curr_step,
        logfile,
        linkcorrlist,
        xyzq,
        pcffile,
        qm_corrdata,
        m2list,
        q1list,
        linkatoms,
        active,
        basedir,
        mminfo,
        qmmminfo,
        pathinfo,
    )
    clean_force = make_clean_force(total_force)

    # function xyz2zmatrix
    # function write Gaussian input with zmat
    # function z2g

    # make new geo
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

    perform_sp(
        gro,
        top,
        xyzq,
        connlist,
        qmatomlist,
        m1list,
        m2list,
        q1list,
        qmmmtop,
        qminfo,
        mminfo,
        qmmminfo,
        linkcorrlist,
        flaglist,
        pcffile,
        linkatoms,
        active,
        logfile,
        basedir,
        step,
        pathinfo,
    )

    return done


def perform_nma(qmmmInputs):

    logger(logfile, "------This will be a numerical) normal mode analysis.------\n")
    logger(
        logfile,
        "Generating a numerical Hessian for the active region using a displacement step of "
        + str(qmmminfo[6])
        + " a.u.\n",
    )
    logger(
        logfile,
        "Will require " + str(len(active) * 6 + 1) + " single point calculations!\n",
    )
    qmenergy, mmenergy, qm_corrdata = get_energy(
        gro,
        qmmminfo[0],
        xyzq,
        connlist,
        qmatomlist,
        m1list,
        m2list,
        q1list,
        qmmmtop,
        qminfo,
        mminfo,
        qmmminfo,
        linkcorrlist,
        flaglist,
        pcffile,
        linkatoms,
        0,
        logfile,
        basedir,
        pathinfo,
    )
    start_energy = float(qmenergy) + float(mmenergy)
    start_forces = read_forces(
        qmatomlist,
        m1list,
        qmmmtop,
        qminfo,
        qmmminfo[0],
        0,
        logfile,
        linkcorrlist,
        xyzq,
        pcffile,
        qm_corrdata,
        m2list,
        q1list,
        linkatoms,
        active,
        basedir,
        mminfo,
        qmmminfo,
        pathinfo,
    )
    start_grad = np.array(start_forces) * -1.0
    hessian_xyz_full = []
    for curr_atom in active:
        grad_deriv_vec = nma_stuff.get_xyz_2nd_deriv(
            curr_atom,
            start_energy,
            start_grad,
            gro,
            top,
            xyzq,
            connlist,
            qmatomlist,
            m1list,
            m2list,
            q1list,
            qmmmtop,
            qminfo,
            mminfo,
            qmmminfo,
            linkcorrlist,
            flaglist,
            pcffile,
            linkatoms,
            active,
            logfile,
            basedir,
        )
        hessian_xyz_full.extend(grad_deriv_vec)
    prep_hess = nma_stuff.prepare_hess(hessian_xyz_full, active)
    evals = []
    for i in range(0, len(prep_hess[0]) - 6):
        evals.extend([float(1000.0)])
    nma_stuff.write_pseudofchk_file(
        qmmminfo[0], evals, prep_hess, prep_hess, active, qmmmtop, logfile, xyzq
    )  # using prep_hess as pseudo-nm_matrix since we do not know yet
    logger(logfile, "Wrote pseudofchk (G03 format).\n")
    write_hess.hes_xyz_fchk(
        str(qmmminfo[0] + ".pseudofchk"), str(qmmminfo[0] + ".hess")
    )
    logger(logfile, "Wrote orca format .hess file.\n")
    evals, nm_matrix = nma.nma_3Nminus6dof_asfunction(
        str(qmmminfo[0] + ".hess"), basedir
    )
    print(nma_stuff.log_nma(
        qmmminfo, logfile, evals, nm_matrix, active, qmmmtop, xyzq, prep_hess
    ))


def run_g16(qmfile, qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    g16cmd = qmmmInputs.pathinfo.g16cmd
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    
    insert = ""
    if int(curr_step) != 0:
        insert = str("." + str(int(curr_step)))


    if not os.path.isfile(str(qmfile) + ".log"):
        logger(logfile, "Running G16 file.\n")
        subprocess.call([g16cmd, str(qmfile)])
        logname = qmfile[:-3]
        logname += "log"
        subprocess.call(["mv", logname, str(jobname + insert + ".gjf.log")])
        subprocess.call(["mv", "fort.7", str(jobname + insert + ".fort.7")])
        logger(logfile, "G16 done.\n")
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
    prefix =  qmmmInputs.pathinfo.gmxpath + qmmmInputs.pathinfo.gmxcmd
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


def get_energy(edrname, qmmmInputs):
    qmenergy, qm_corrdata = get_qmenergy(qmfile, qmmmInputs)
    mmenergy = get_mmenergy(str(edrname), qmmmInputs)
    linkcorrenergy = get_linkenergy_au(qm_corrdata, qmmmInputs)

    qmenergy -= linkcorrenergy
    methodstring = str(qmmmInputs.qmparams.method)
    if qminfo[2] != "NONE":
        methodstring += str("/" + str(qmmmInputs.qmparams.basis))
    logger(
        logfile,
        str(
            "Single point energy done. QM/MM energy is {:>20.10f} (QM, link atom corrected ".format(
                float(qmenergy)
            )
            + methodstring
            + ") + {:>20.10f} (MM) = {:>20.10f} (a.u.)\n".format(
                float(mmenergy), float(qmenergy) + float(mmenergy)
            )
        ),
    )

    return qmenergy, mmenergy, linkcorrenergy, qm_corrdata

def get_qmenergy(qmfile, qmmmInputs, pcffile):
    qmprog = qmmmInputs.qmparams.program
    extra_string = qmmmInputs.qmparams.extra
    pcffile = qmmmInputs.qmmmparams.pcffile
    logfile = qmmmInputs.logfile
    basedir = qmmmInputs.qmmmparams.basedir

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
                    logger(logfile, "Obtaining charge self-interaction...")
                    pcf_self_pot = read_pcf_self(qmfile)
                    logger(
                        logfile, "done: {:>20.10f} a.u.\n".format(float(pcf_self_pot))
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
    prefix =  qmmmInputs.pathinfo.gmxpath + qmmmInputs.pathinfo.gmxcmd
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
    xyzq = qmmmInput.qmmmparams.xyzq
    linkcorrlist = qmmmInput.qmmmparams.linkcorrlist
    m1list = qmmmInput.qmmmparams.m1list
    m2list = qmmmInput.qmmmparams.m2list
    q1list = qmmmInput.qmmmparams.q1list
    qmmmtop = qmmmInput.qmmmparams.qmmmtop
    linkatoms = qmmmInput.qmmmparams.linkatoms
    logfile = qmmmInput.qmmmparams.logfile



    pcffile = qmmmInput.qmmmparams.pcffile
    qmatomlist = qmmmInput.qmmmparams.qmatomlist

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
    energycorr = databasecorrection(
        "ENERGY", "aminoacid_CACB", dist, qmmmInputs, basedir, logfile
    )
    linkenergy -= energycorr
    # sign inverted due to correction convention (subtracting)
    return linkenergy



if __name__ == "__main__":
    print("This file serves as a library for gmx2qmmm-related functions.")
    print("Do not execute directly. Use gmx2qmmm instead.")
