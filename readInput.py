from gmx2qmmm.pointcharges import generate_pcf_from_top as make_pcf
from gmx2qmmm.pointcharges import prepare_pcf_for_shift as prep_pcf
from gmx2qmmm.pointcharges import generate_charge_shift as final_pcf
from gmx2qmmm.operations import generate_top as topprep
from gmx2qmmm.operations import qmmm

class QMParams:
    def __init__(self, inp):
        self.inp = inp
        self.program, self.method, self.basis, self.charge, self.multi, self.cores, self.memory, self.extra = self.read_qmparams(inp)
    
    def read_qmparams(self, inp):
        import re
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
        import re
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
        self.g16path, self.tmpath, self.orcapath, self.gmxpath, self.g16cmd, self.tmcmd, self.orcacmd, self.gmxcmd, self.gmxtop = self.read_pathparams(inp)

    def read_pathparams(self, inp):
        import re
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
        self.jobname, self.jobtype, self.propagater, self.maxcycle, self.initstep, self.f_thresh, self.disp, self.curr_step, self.databasefit, self.optlastonly = self.read_qmmmparams(inp)

    def read_qmmmparams(self, inp):
        import re
        jobname="testjob"
        jobtype="singlepoint"
        propa="steep"
        info=[jobname,jobtype.upper(),float(0.00001),int(5),float(0.1),propa.upper(),float(0.0018897261),int(0),"MORSE","YES"]
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
    def __init__(self, inputFiles):
        self.inputFiles = inputFiles
        self.mmparams = MMParams(inputFiles.mmFile)
        self.qmparams = QMParams(inputFiles.qmFile)
        self.pathparams = PathParams(inputFiles.pathFile)
        self.qmmmparams = QMMMParams(inputFiles.qmmmFile)
        self.logfile = inputFiles.logfile


    
        