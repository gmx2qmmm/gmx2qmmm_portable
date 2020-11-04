import re
import os
import numpy as np

from gmx2qmmm._helper import _flatten, logger, stepper
from gmx2qmmm._helper import get_linkatoms_ang, make_xyzq
from gmx2qmmm.pointcharges import generate_pcf_from_top as make_pcf
from gmx2qmmm.pointcharges import prepare_pcf_for_shift as prep_pcf
from gmx2qmmm.pointcharges import generate_charge_shift as final_pcf
from gmx2qmmm.operations import generate_top as topprep
from gmx2qmmm.operations import qmmm

class QMParams:
    def __init__(self, inp):
        self.inp = inp
        (self.program, self.method, self.basis, self.charge, 
            self.multi, self.cores, self.memory, self.extra) = self.read_qmparams(inp)
    
    def read_qmparams(self, inp):
        info = [
            "G16",
            "BP86",
            "STO-3G",
            int(0),
            int(1),
            int(1),
            int(1000),
            "NONE",
        ]  # program method basis charge multiplicity cores memory(MB) extraopts
        with open(inp) as ifile:
            for line in ifile:
                match = re.search(r"^program\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[0] = str(match.group(1)).upper()
                match = re.search(r"^method\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[1] = str(match.group(1)).upper()
                match = re.search(r"^basis\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[2] = str(match.group(1)).upper()
                match = re.search(r"^charge\s*=\s*([-]*\d+)", line, flags=re.MULTILINE)
                if match:
                    info[3] = int(match.group(1))
                match = re.search(r"^multiplicity\s*=\s*(\d+)", line, flags=re.MULTILINE)
                if match:
                    info[4] = int(match.group(1))
                match = re.search(r"^cores\s*=\s*(\d+)", line, flags=re.MULTILINE)
                if match:
                    info[5] = int(match.group(1))
                match = re.search(r"^memory\s*=\s*(\d+)", line, flags=re.MULTILINE)
                if match:
                    info[6] = int(match.group(1))
                match = re.search(r"^extra\s*=\s*(.+)$", line, flags=re.MULTILINE)
                if match:
                    info[7] = match.group(1)
        return info

class MMParams:
    def __init__(self, inp):
        self.inp = inp
        info, flaglist = self.read_mmparams(inp)
        self.fField, self.rvdw = info
        self.flaglist = flaglist
        
    def read_mmparams(self, inp):
        info = ["amberGS", float(2.0)]
        flaglist = []
        with open(inp) as ifile:
            for line in ifile:
                match = re.search(r"^ff\s*=\s*(\d*\.*\d*)", line, flags=re.MULTILINE)
                if match:
                    info[0] = float(match.group(1))
                    continue
                match = re.search(r"^rvdw\s*=\s*(\d*\.*\d*)", line, flags=re.MULTILINE)
                if match:
                    info[1] = float(match.group(1))
                    continue
                match = re.search(r"\-D(\S*)", line, flags=re.MULTILINE)
                if match:
                    if match.group(1) not in flaglist:
                        flaglist.append((match.group(1)).upper())

        return info, flaglist

class PathParams:
    def __init__(self, inp):
        self.inp = inp
        (self.g16path, self.tmpath, self.orcapath, self.gmxpath,
         self.g16cmd, self.tmcmd, self.orcacmd, self.gmxcmd,
         self.gmxtop) = self.read_pathparams(inp)

    def read_pathparams(self, inp):
        # software path
        g16path = ""
        tmpath = ""
        orcapath = ""
        gmxpath = "gromacs-2019.1/bin/"
        # Gaussain,Turbomole,ORCA, GROMACS, gaussianCMD, TMCMD, orcaCMD, gmxCMD, gmxTop
        info = [
            "",
            "",
            "",
            "gromacs-2019.1/bin/",
            "rung16",
            "",
            "",
            "gmx19",
            "gromacs-2019.1/share/gromacs/top/",
        ]
        with open(inp) as ifile:
            for line in ifile:
                match = re.search(r"^g16path\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[0] = str(match.group(1))

                match = re.search(r"^tmpath\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[1] = str(match.group(1))

                match = re.search(r"^orcapath\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[2] = str(match.group(1))

                match = re.search(r"^gmxpath\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[3] = str(match.group(1))

                match = re.search(r"^g16cmd\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[4] = str(match.group(1))

                match = re.search(r"^tmcmd\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[5] = str(match.group(1))

                match = re.search(r"^orcacmd\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[6] = str(match.group(1))

                match = re.search(r"^gmxcmd\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[7] = str(match.group(1))

                match = re.search(r"^gmxtop_path\s*=\s*(\S+)", line, flags=re.MULTILINE)
                if match:
                    info[8] = str(match.group(1))
        return info

class QMMMParams:
    def __init__(self, inp):
        self.inp = inp
        (self.jobname, self.jobtype, self.propagater, self.maxcycle, self.initstep, self.f_thresh, 
            self.disp, self.curr_step, self.databasefit, self.optlastonly) = self.read_qmmmparams(inp)

    def read_qmmmparams(self, inp):
        jobname="testjob"
        jobtype="singlepoint"
        propa="steep"
        info=[jobname,jobtype.upper(),float(0.00001),int(5),float(0.1),
                propa.upper(),float(0.0018897261),int(0),"MORSE","YES"]
        with open(inp) as ifile:
            for line in ifile:
                match=re.search(r'^jobname\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[0]=match.group(1)
                match=re.search(r'^jobtype\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[1]=str(match.group(1)).upper()
                match=re.search(r'^f_thresh\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[2]=float(match.group(1))
                match=re.search(r'^maxcycle\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[3]=int(match.group(1))
                match=re.search(r'^initstep\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[4]=float(match.group(1))
                match=re.search(r'^propagator\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[5]=str(match.group(1)).upper()
                match=re.search(r'^disp\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[6]=float(match.group(1))
                match=re.search(r'^current_step\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[7]=int(match.group(1))
                match=re.search(r'^databasefit\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    #info[8]=match.group(1)
                    info[8]=str(match.group(1)).upper()
                match=re.search(r'^optlastonly\s*=\s*(\S+)', line,flags=re.MULTILINE)
                if match:
                    info[9]=str(match.group(1)).upper()
        return info

class QMMMInputs:
    """docstring for ClassName"""
    def __init__(self, inputFiles, basedir):
        self.inputFiles = inputFiles
        self.basedir = basedir
        self.logfile = inputFiles.logfile
        logfile = self.logfile
        #Read .dat file
        logger(logfile, "Reading input parameters...\n")
        self.mmparams = MMParams(inputFiles.mmFile)
        self.qmparams = QMParams(inputFiles.qmFile)
        self.pathparams = PathParams(inputFiles.pathFile)
        self.qmmmparams = QMMMParams(inputFiles.qmmmFile)
        logger(logfile, "Done.\n")

        logger(logfile, "Initializing dependencies...")
        self.gro = stepper(inputFiles.coord, self.qmmmparams.curr_step)
        logger(self.logfile, "complete.\n")


        self.chargevec = []
        logger(logfile, "Trying to understand your MM files.\n")
        logger(logfile, "List of molecules...")
        self.mollist = make_pcf.readmols(inputFiles.top)
        logger(logfile, "done.\n")
        logger(logfile, "Reading charges...")
        for element in self.mollist:
            self.chargevec.extend(make_pcf.readcharges(element, inputFiles.top, self.pathparams.gmxtop))
        logger(logfile, "done.\n")

        
        if inputFiles.coord[-4:] == ".gro":
            logger(logfile, "Reading geometry (.gro)...")
            self.geo = make_pcf.readgeo(inputFiles.coord)
        elif inputFiles.coord[-4:] == ".g96":
            logger(logfile, "Reading geometry (.g96)...")
            self.geo = make_pcf.readg96(inputFiles.coord)
        logger(logfile, "%s\n" % self.geo)
        logger(logfile, "done.\n")

        logger(logfile, "Reading connectivity matrix...")
        self.connlist = read_conn_list_from_top(inputFiles.top, self.mollist, self.basedir)
        logger(logfile, "done.\n")
        
        logger(logfile, "Trying to understand your QM, MM and QM/MM parameters.\n")
        logger(logfile, "Reading QM atom list...")
        self.qmatomlist = prep_pcf.read_qmatom_list(inputFiles.qmatoms)
        logger(logfile, "done.\n")

        if self.gro[-4:] == ".gro":
            logger(logfile, "Writing high-precision coordinate file...")
            grohigh = write_highprec(self.gro, self.qmmmparams.jobname, logfile)
            self.gro = self.qmmmparams.jobname + ".g96"
            logger(logfile, "done.\n")

        logger(
            logfile,
            'Setting up a calculation named "'
            + str(self.qmmmparams.jobname)
            + '" of type "'
            + str(self.qmmmparams.jobtype)
            + '".\n',
        )

        self.qmmmtop = str(self.qmmmparams.jobname + ".qmmm.top")
        self.pcffile = str(self.qmmmparams.jobname + ".pointcharges")
        self.pcffile = stepper(self.pcffile, self.qmmmparams.curr_step)

        logger(
            logfile,
            "Starting the preparation of the point charge field used in the calculation.\n",
        )
        logger(logfile, "Creating the full xyzq Matrix...")
        self.xyzq = make_xyzq(self.geo, self.chargevec)
        logger(logfile, "done.\n")
        
        logger(
            logfile,
            "Preparing the point charge field for a numerically optimized charge shift...",
        )
        (
            self.qmcoordlist,
            self.m1list,
            self.m2list,
            self.updated_chargelist,
        ) = prep_pcf.prepare_pcf_for_shift_fieldsonly(self.xyzq, self.qmatomlist, self.qmmmparams.maxcycle, self.connlist)
        logger(logfile, "done.\n")

        logger(logfile, "Setting up link atoms...")
        self.linkatoms = get_linkatoms_ang(self.xyzq, self.qmatomlist, self.m1list, self.connlist, [])
        self.linkcorrlist, self.q1list, self.q2list, self.q3list, self.m3list = get_linkcorrlist(
            self.linkatoms, self.qmatomlist, self.m1list, self.m2list, self.connlist
        )
        logger(logfile, "done.\n")

        if self.qmmmparams.curr_step > 0:
            self.qmmmparams.jobname = str(self.qmmmparams.jobname + "." + str(self.qmmmparams.curr_step))
        if not os.path.isfile(str(self.qmmmparams.jobname + ".pointcharges")):
            logger(logfile, "Shifting...")
            final_pcf.generate_charge_shift_fieldsonly(
                self.updated_chargelist, self.m1list, self.qmcoordlist, self.m2list, self.qmmmparams.jobname, self.basedir
            )
        else:
            logger(
                logfile,
                "NOTE: Shifting omitted due to "
                + str(self.qmmmparams.jobname + ".pointcharges")
                + " being an existing file!\n",
            )
        logger(logfile, "done.\n")

        logger(logfile, "Preparing the QM/MM top file...")
        topprep.generate_top_listsonly(
            self.inputFiles.top,
            self.inputFiles.qmatoms,
            self.qmmmtop,
            self.mmparams.flaglist,
            self.q1list,
            self.q2list,
            self.q3list,
            self.m1list,
            self.m2list,
            self.m3list,
            self.basedir,
            self.logfile,
            self.pathparams.gmxtop,
        )
        logger(logfile, "done.\n")
        active = []
        if jobtype != "SINGLEPOINT":
            logger(logfile, "Reading indices of active atoms...")
            active = prep_pcf.read_qmatom_list(inputFiles.act)
            logger(logfile, "done.\n")
 
def get_curr_top(molname, top, basedir):
    curr_top = top
    found = make_pcf.checkformol(molname, top)

    if not found:
        toplist = make_pcf.getincludelist(top)
        for element in toplist:
            found = make_pcf.checkformol(molname, element)
            if found:
                curr_top = element
                break
    if not found:
        print("No charges found for " + str(molname) + ". Exiting.")
        exit(1)
    return curr_top

def get_mollength_direct(molname, top):
    mollength = 0
    with open(top) as ifile:
        # print str(top) + " is open"
        for line in ifile:
            match = re.search(r"^;", line, flags=re.MULTILINE)
            if match:
                continue
            match = re.search(r"^\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
            if match:
                # print "moltype was found"
                for line in ifile:
                    match = re.search(r"^;", line, flags=re.MULTILINE)
                    if match:
                        continue
                    else:
                        matchstring = r"^\s*" + re.escape(molname)
                        match = re.search(matchstring, line, flags=re.MULTILINE)
                        if match:
                            # print str(matchstring) + " was found"
                            for line in ifile:
                                match = re.search(r"^;", line, flags=re.MULTILINE)
                                if match:
                                    continue
                                else:
                                    match = re.search(
                                        r"^\[\s*atoms\s*\]", line, flags=re.MULTILINE
                                    )
                                    if match:
                                        # print "atoms was found"
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

def get_connlist(offset, molname, top):
    connlist = []
    with open(top) as ifile:
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
                        matchstring = r"^\s*" + re.escape(molname)
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
                                                    for element in connlist:
                                                        if int(element[0]) == a:
                                                            if int(b) not in np.array(
                                                                element
                                                            ).astype(int):
                                                                element.append(int(b))
                                                                found = True
                                                    if not found:
                                                        connlist.append(
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
    return connlist

def read_conn_list_from_top(top, mollist, basedir):
    count = 0
    connlist = []
    for element in mollist:
        curr_top = get_curr_top(element[0], top, basedir)
        for i in range(0, int(element[1])):
            mollength = get_mollength_direct(element[0], curr_top)
            connset = get_connlist(count, element[0], curr_top)
            connlist += connset
            count += int(mollength)
    return connlist

def get_linkcorrlist(linkatoms, qmatomlist, m1list, m2list, connlist):
    linkcorrlist = []
    m3list = []
    m4list = []
    q1list = []
    q2list = []
    q3list = []
    # get q1 and m2
    for element in m1list:
        q1line = []
        for entry in connlist:
            if int(element) in np.array(entry).astype(int):
                if int(element) != int(entry[0]):
                    if int(entry[0]) in np.array(qmatomlist).astype(int):
                        q1line.append(int(entry[0]))
                else:
                    for i in range(1, len(entry)):
                        if int(entry[i]) in np.array(qmatomlist).astype(int):
                            q1line.append(int(entry[i]))
        q1list.append(q1line)
    # get q2
    q1list = list(_flatten(q1list))
    for element in q1list:
        q2line = []
        for conn in connlist:
            if int(element) in np.array(conn).astype(int):
                if (
                    int(element) != int(conn[0])
                    and (int(conn[0]) in np.array(qmatomlist).astype(int))
                    and (int(conn[0]) not in np.array([int(x) for x in _flatten(q1list)]))
                ):
                    q2line.append(int(conn[0]))
                elif int(element) == int(conn[0]):
                    for i in range(1, len(conn)):
                        if (int(conn[i]) in np.array(qmatomlist).astype(int)) and (
                            int(conn[i]) not in np.array([int(x) for x in _flatten(q1list)])
                        ):
                            q2line.append(int(conn[i]))
        q2list.append(q2line)
    # get q3
    for element in q2list:
        q3lineline = []
        for entry in element:
            q3line = []
            for conn in connlist:
                if int(entry) in np.array(conn).astype(int):
                    if (
                        int(entry) != int(conn[0])
                        and (int(conn[0]) in np.array(qmatomlist).astype(int))
                        and int(conn[0]) not in np.array(q1list).astype(int)
                        and int(conn[0]) not in np.array([int(x) for x in _flatten(q1list)])
                    ):
                        q3line.append(int(conn[0]))
                    elif int(entry) == int(conn[0]):
                        for i in range(1, len(conn)):
                            if (
                                int(conn[i]) in np.array(qmatomlist).astype(int)
                                and int(conn[i]) not in np.array(q1list).astype(int)
                                and int(conn[i]) not in np.array([int(x) for x in _flatten(q1list)])
                            ):
                                q3line.append(int(conn[i]))
            q3lineline.append(q3line)
        q3list.append(q3lineline)
    # get m3
    for element in m2list:
        m3lineline = []
        for entry in element:
            m3line = []
            for conn in connlist:
                if int(entry) in np.array(conn).astype(int):
                    if (
                        int(entry) != int(conn[0])
                        and (int(conn[0]) not in np.array(qmatomlist).astype(int))
                        and int(conn[0]) not in np.array(m1list).astype(int)
                    ):
                        m3line.append(int(conn[0]))
                    elif int(entry) == int(conn[0]):
                        for i in range(1, len(conn)):
                            if int(conn[i]) not in np.array(qmatomlist).astype(
                                int
                            ) and int(conn[i]) not in np.array(m1list).astype(int):
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
                for conn in connlist:
                    if int(stuff) in np.array(conn).astype(int):
                        if (
                            int(stuff) != int(conn[0])
                            and (int(conn[0]) not in np.array(qmatomlist).astype(int))
                            and int(conn[0]) not in np.array(m1list).astype(int)
                        ):
                            found = False
                            for morestuff in m2list:
                                if int(conn[0]) in np.array(morestuff).astype(int):
                                    found = True
                                    break
                            if not found:
                                m4line.append(int(conn[0]))
                        elif int(stuff) == int(conn[0]):
                            for i in range(1, len(conn)):
                                if int(conn[i]) not in np.array(qmatomlist).astype(
                                    int
                                ) and int(conn[i]) not in np.array(m1list).astype(int):
                                    found = False
                                    for morestuff in m2list:
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
    count = 0
    for element in m1list:
        linkpairline = []
        for entry in m2list[count]:
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
        for stuff in m2list[count]:
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
            for stuff in m2list[count]:
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