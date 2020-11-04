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
    qmprogram = qmmmInputs.qmparams
    logfile = qmmmInputs.logfile

    if qmprogram == "G16":
        make_g16_inp(qmmmInputs)
        qmenergy, mmenergy, qm_corrdata = get_energy(qmmmInputs)
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
    oforce = open("oforce.txt", "w")
    for i in range(len(total_force)):
        oforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
        )
    oforce.close()

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


if __name__ == "__main__":
    print("This file serves as a library for gmx2qmmm-related functions.")
    print("Do not execute directly. Use gmx2qmmm instead.")
