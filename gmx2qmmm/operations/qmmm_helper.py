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

def qmmm_prep(qmmmInputs, new_gro, top,qminfo,qmatomlist, connlist, linkatoms,  basedir,  logfile, pathinfo):
    jobname = qmmmInputs.qmmmparams.jobname
    curr_step = qmmmInputs.qmmmparams.curr_step
    top = qmmmInputs.top
    qmatomlist = qmmmInputs.qmatomlist
    linkatoms = qmmmInputs.linkatoms
    basedir = qmmmInputs.basedir
    logfile = qmmmInputs.logfile
    charge = qmmmInputs.qmparams.charge
    gmxtop_path = qmmmInputs.pathparams.gmxtops

    geo = make_pcf.readg96(new_gro)
    logger(logfile, "List of molecules...")
    mollist = make_pcf.readmols(top)
    logger(logfile, "done.\n")
    logger(logfile, "Reading charges...")
    chargevec = []
    for element in mollist:
        chargevec.extend(make_pcf.readcharges(element, top, gmxtop_path))
    logger(logfile, "done.\n")
    new_xyzq = make_xyzq(geo, chargevec)
    logger(logfile, str("Made new xyzq matrix.\n"))
    logger(
        logfile,
        "Preparing the point charge field for a numerically optimized charge shift...",
    )
    (
        qmcoordlist,
        m1list,
        m2list,
        updated_chargelist,
    ) = prep_pcf.prepare_pcf_for_shift_fieldsonly(
        new_xyzq, qmatomlist, charge, connlist
    )
    logger(logfile, "done.\n")
    new_links = get_linkatoms_ang(
        new_xyzq, qmatomlist, m1list, connlist, linkatoms
    )
    logger(logfile, str("Updated positions of link atoms.\n"))
    filename = jobname
    if curr_step != 0:
        filename += "." + str(int(curr_step))
    if not os.path.isfile(str(filename + ".pointcharges")):
        logger(logfile, "Shifting...")
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
    logger(logfile, "done.\n")
    return new_xyzq, m1list, m2list, new_links


if __name__ == "__main__":
    print("This file serves as a library for gmx2qmmm-related functions.")
    print("Do not execute directly. Use gmx2qmmm instead.")
