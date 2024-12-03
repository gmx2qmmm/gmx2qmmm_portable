#   // INITIAL DESCRIPTION //
"""Short Module Description; Reference To Readme"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-09-25'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import sqlite3
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from gmx2qmmm.logging import Logger
from gmx2qmmm.generators._helper import filter_xyzq, _flatten


#   // TODOS & NOTES //
#   TODO:
#   - change 'make_xyzq' and 'make_xyzq_io' to be one function
#   NOTE: For now I'm just randomly keeping functions regarding geometry stuff here. Sorting into a class later (AJ)

class GeneratorQMMM():

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, work_dir, base_dir, class_qm_job, class_mm_job) -> None:
        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.work_dir = work_dir
        self.base_dir = base_dir
        self.class_qm_job = class_qm_job
        self.class_mm_job = class_mm_job

    def get_m2charges(self):
        # XX AJ check if this function is necessary
        # XX AJ because of lazyness we have this function in each class in this file! Change that, Alina!
        self.m2charges = []
        count = 0
        for _ in self.system.list_atoms_m1:
            m2chargeline = []
            for i in range(0, len(self.system.list_atoms_m2[count])):
                m2chargeline.append(float(self.system.array_xyzq_current[int(self.system.list_atoms_m2[count][i]) - 1][3]))
            count += 1
            self.m2charges.append(m2chargeline)

    def read_pcffile(self):
            # XX AJ because of lazyness we have this function in each class in this file! Change that, Alina!

            '''
            ------------------------------ \\
            EFFECT: \\
            --------------- \\
            XX \\
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

            pcf = []
            with open(self.PCF.pcf_filename) as ifile:
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

    def databasecorrection(self, energy_or_force, cut, distance):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        XX \\
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

        # XX AJ check what this connection object is
        # XX AJ: I had to convert the path object to string to get the sqlite object, is there a better way to do that?
        self.base_dir = str(self.base_dir)
        conn = sqlite3.connect(self.base_dir + "/correction_database/database.sqlite")

        # check if method exist in database
        method_set = conn.cursor()
        method_set.execute(
            "SELECT * FROM "
            + cut
            + ' WHERE forcefield="'
            + self.dict_input_userparameters['forcefield']
            + '" AND method="'
            + self.dict_input_userparameters['method']
            + '" AND basisset ="'
            + self.dict_input_userparameters['basis']
            + '"'
        )
        method_set_value = method_set.fetchall()
        if len(method_set_value) == 0:
            cut = "aminoacid_CACB"
            self.dict_input_userparameters['forcefield'] = "amberGS"
            self.dict_input_userparameters['method'] = "CAM-B3LYP"
            self.dict_input_userparameters['basis'] = "6-31++G**"
            # logger(
            #     logfile,
            #     "Unexisted method in correction database, changing to default correction method...\n",
            # )

        c = conn.cursor()
        c.execute(
            "SELECT * FROM "
            + cut
            + ' WHERE forcefield="'
            + self.dict_input_userparameters['forcefield']
            + '" AND method="'
            + self.dict_input_userparameters['method']
            + '" AND basisset ="'
            + self.dict_input_userparameters['basis']
            + '"'
        )
        db_values = c.fetchall()[0]

        conn.close()
        returnvalue = 0
        if len(db_values) > 0:
            if energy_or_force == "energy":
                if self.dict_input_userparameters['databasefit'] == "poly":
                    returnvalue = (
                        db_values[5] * distance * distance * distance
                        + db_values[6] * distance * distance
                        + db_values[7] * distance
                        + db_values[8]
                    )
                elif self.dict_input_userparameters['databasefit'] == "morse":
                    returnvalue = (
                        db_values[9]
                        * (
                            np.exp(-2 * db_values[10] * (distance - db_values[11]))
                            - 2 * np.exp(-db_values[10] * (distance - db_values[11]))
                        )
                        + db_values[12]
                    )
                elif self.dict_input_userparameters['databasefit'] == "no":
                    returnvalue = 0
                    # logger(logfile, "No energy correction.\n")
            elif energy_or_force == "forces":
                if self.dict_input_userparameters['databasefit'] == "poly":
                    returnvalue = (
                        db_values[13] * distance * distance * distance
                        + db_values[14] * distance * distance
                        + db_values[15] * distance
                        + db_values[16]
                    )
                elif self.dict_input_userparameters['databasefit'] == "morse":
                    returnvalue = (
                        db_values[17]
                        * (
                            np.exp(-2 * db_values[18] * (distance - db_values[19]))
                            - 2 * np.exp(-db_values[18] * (distance - db_values[19]))
                        )
                        + db_values[20]
                    )
                elif self.dict_input_userparameters['databasefit'] == "no":
                    returnvalue = 0
                    # logger(logfile, "No force correction.\n")
        return returnvalue


class GeneratorEnergies(GeneratorQMMM):

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, work_dir, base_dir, class_qm_job, class_mm_job) -> None:
        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.work_dir = work_dir
        self.base_dir = base_dir
        self.class_qm_job = class_qm_job
        self.class_mm_job = class_mm_job

    def get_linkcorrenergy(self):
        if self.system.linkcorrlist:
            linkcorrenergy = self.get_linkenergy_au()
        else:
            linkcorrenergy = 0.0
        return linkcorrenergy


    def get_linkenergy_au(self):
        linkenergy = 0.0
        self.get_m2charges() # XX AJ correct
        for element in self.system.linkcorrlist:
            z1 = 0.0
            v1 = []
            v2 = []
            if int(element[0]) in np.array(list(_flatten(self.system.list_atoms_m2))).astype(int):
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        if int(self.system.list_atoms_m2[i][j]) == int(element[0]):
                            z1 = float(self.m2charges[i][j])
                            v1 = [
                                self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                                self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                                self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                            ]
                            break
                    if z1 != 0.0:
                        break
            elif int(element[0]) in np.array(self.system.list_atoms_qm).astype(int):
                for i in range(0, len(self.system.list_atoms_qm)):
                    if int(self.system.list_atoms_qm[i]) == int(element[0]):
                        z1 = float(self.class_qm_job.qm_corrdata[i][2])
                        v1 = [
                            self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                            self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                            self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                        ]
                        break
            elif int(element[0]) in np.array(self.system.list_atoms_m1).astype(int):
                for i in range(0, len(self.system.list_atoms_m1)):
                    if int(self.system.list_atoms_m1[i]) == int(element[0]):
                        z1 = float(self.class_qm_job.qm_corrdata[i + len(self.system.list_atoms_qm)][2])
                        v1 = [
                            self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
                        ]
                        break
            else:
                z1 = float(self.system.array_xyzq_current[int(element[0]) - 1][3])
                v1 = [
                    self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                    self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                    self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                ]
            z2 = 0.0
            if int(element[1]) in _flatten(self.system.list_atoms_m2):
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        if int(self.system.list_atoms_m2[i][j]) == int(element[1]):
                            z2 = float(self.m2charges[i][j])
                            v2 = [
                                self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                                self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                                self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                            ]
                            break
                    if z2 != 0.0:
                        break
            elif int(element[1]) in np.array(self.system.list_atoms_qm).astype(int):
                for i in range(0, len(self.system.list_atoms_qm)):
                    if int(self.system.list_atoms_qm[i]) == int(element[1]):
                        z2 = float(self.class_qm_job.qm_corrdata[i][2])
                        v2 = [
                            self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                            self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                            self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                        ]
                        break
            elif int(element[1]) in np.array(self.system.list_atoms_m1).astype(int):
                for i in range(0, len(self.system.list_atoms_m1)):
                    if int(self.system.list_atoms_m1[i]) == int(element[1]):
                        z2 = float(self.class_qm_job.qm_corrdata[i + len(self.system.list_atoms_qm)][2])
                        v2 = [
                            self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
                        ]
                        break
            else:
                z2 = float(self.system.array_xyzq_current[int(element[1]) - 1][3])
                v2 = [
                    self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                    self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                    self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                ]
            v12 = np.array(v1) - np.array(v2)
            dist = np.linalg.norm(v12)
            linkenergy += z1 * z2 / dist
        # now also all atoms in the corrdata list with the mod and linkcorr point charges
        # mod first. mod is charge in pcffile minus m2charge
        pcf = self.read_pcffile()
        for i in range(0, len(self.system.list_atoms_m2)):
            for j in range(0, len(self.system.list_atoms_m2[i])):
                curr_mod = []
                for k in range(0, 3):
                    curr_mod.append(float(pcf[int(self.system.list_atoms_m2[i][j]) - 1][k]) / 0.52917721)
                curr_mod_charge = (
                    float(float(pcf[int(self.system.list_atoms_m2[i][j]) - 1][3])) - self.m2charges[i][j]
                )
                for k in range(0, len(self.system.list_atoms_qm)):
                    v1 = [
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][0] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][1] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][2] / 0.52917721,
                    ]
                    z1 = float(self.class_qm_job.qm_corrdata[k][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    linkenergy += z1 * curr_mod_charge / dist
                for k in range(0, len(self.system.list_coordinates_linkatoms)):
                    v1 = [
                        self.system.list_coordinates_linkatoms[k][0] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][1] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][2] / 0.52917721,
                    ]
                    z1 = float(self.class_qm_job.qm_corrdata[k + len(self.system.list_atoms_qm)][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    linkenergy += z1 * curr_mod_charge / dist
        # now linkcorr. linkcorr are last m2*2 entries in pcf
        m2count = 0
        linkstart = len(pcf) - 2 * len(list(_flatten(self.system.list_atoms_m2)))
        for i in range(0, len(self.system.list_atoms_m2)):
            for j in range(0, len(self.system.list_atoms_m2[i])):
                curr_mod = []
                for k in range(0, 3):
                    curr_mod.append(float(pcf[int(linkstart) + m2count][k]) / 0.52917721)
                curr_mod_charge = float(float(pcf[int(linkstart) + m2count][3]))
                m2count += 1
                for k in range(0, len(self.system.list_atoms_qm)):
                    v1 = [
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][0] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][1] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][2] / 0.52917721,
                    ]
                    z1 = float(self.class_qm_job.qm_corrdata[k][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    linkenergy += z1 * curr_mod_charge / dist
                for k in range(0, len(self.system.list_coordinates_linkatoms)):
                    v1 = [
                        self.system.list_coordinates_linkatoms[k][0] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][1] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][2] / 0.52917721,
                    ]
                    z1 = float(self.class_qm_job.qm_corrdata[k + len(self.system.list_atoms_qm)][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    linkenergy += z1 * curr_mod_charge / dist
        # now, add the correction of energy for the link atoms. currently only C-C bond cuts supported.
        for i in range(0, len(self.system.list_coordinates_linkatoms)):
            v1 = [
                self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
            ]
            _flattened = list(_flatten(self.system.list_atoms_q1))
            v2 = [
                self.system.array_xyzq_current[int(_flattened[i]) - 1][0] / 0.52917721,
                self.system.array_xyzq_current[int(_flattened[i]) - 1][1] / 0.52917721,
                self.system.array_xyzq_current[int(_flattened[i]) - 1][2] / 0.52917721,
            ]
            v12 = np.array(v2) - np.array(v1)
            dist = np.linalg.norm(v12)
        dist = dist * 0.7409471631
        energycorr = self.databasecorrection("energy", "aminoacid_CACB", dist)
        linkenergy -= energycorr
        # sign inverted due to correction convention (subtracting)
        return linkenergy

class GeneratorForces(GeneratorQMMM):

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, work_dir, base_dir, class_qm_job, class_mm_job) -> None:
        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.work_dir = work_dir
        self.base_dir = base_dir
        self.class_qm_job = class_qm_job
        self.class_mm_job = class_mm_job

    def get_linkcorrforces(self):
        if self.system.linkcorrlist:
            linkcorrforces = self.get_linkforces_au()
        else:
            linkcorrforces = 0.0
        return linkcorrforces

    def get_linkforces_au(self):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        XX \\
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

        linkforces = []
        # Force Coulomb: z1*z2*(distance along coord)/(distance between charges)**3
        for element in self.system.array_xyzq_current:  # this is just to count an entry for each atom!
            linkforces.append([0.0, 0.0, 0.0])
        m2charges = self.get_m2charges()
        for element in self.system.linkcorrlist:
            z1 = 0.0
            v1 = []
            v2 = []
            if int(element[0]) in _flatten(self.system.list_atoms_m2):
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        if int(self.system.list_atoms_m2[i][j]) == int(element[0]):
                            z1 = float(m2charges[i][j])
                            v1 = [
                                self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                                self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                                self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                            ]
                            break
                    if z1 != 0.0:
                        break
            elif int(element[0]) in np.array(self.system.list_atoms_qm).astype(int):
                for i in range(0, len(self.system.list_atoms_qm)):
                    if int(self.system.list_atoms_qm[i]) == int(element[0]):
                        z1 = float(self.class_qm_job.qm_corrdata[i][2])
                        v1 = [
                            self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                            self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                            self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                        ]
                        break
            elif int(element[0]) in np.array(self.system.list_atoms_m1).astype(int):
                for i in range(0, len(self.system.list_atoms_m1)):
                    if int(self.system.list_atoms_m1[i]) == int(element[0]):
                        z1 = float(self.class_qm_job.qm_corrdata[i + len(self.system.list_atoms_qm)][2])
                        v1 = [
                            self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
                        ]
                        break
            else:
                z1 = float(self.system.array_xyzq_current[int(element[0]) - 1][3])
                v1 = [
                    self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                    self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                    self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                ]
            z2 = 0.0
            if int(element[1]) in _flatten(self.system.list_atoms_m2):
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        if int(self.system.list_atoms_m2[i][j]) == int(element[1]):
                            z2 = float(self.m2charges[i][j])
                            v2 = [
                                self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                                self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                                self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                            ]
                            break
                    if z2 != 0.0:
                        break
            elif int(element[1]) in np.array(self.system.list_atoms_qm).astype(int):
                for i in range(0, len(self.system.list_atoms_qm)):
                    if int(self.system.list_atoms_qm[i]) == int(element[1]):
                        z2 = float(self.class_qm_job.qm_corrdata[i][2])
                        v2 = [
                            self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                            self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                            self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                        ]
                        break
            elif int(element[1]) in np.array(self.system.list_atoms_m1).astype(int):
                for i in range(0, len(self.system.list_atoms_m1)):
                    if int(self.system.list_atoms_m1[i]) == int(element[1]):
                        z2 = float(self.class_qm_job.qm_corrdata[i + len(self.system.list_atoms_qm)][2])
                        v2 = [
                            self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
                        ]
                        break
            else:
                z2 = float(self.system.array_xyzq_current[int(element[1]) - 1][3])
                v2 = [
                    self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                    self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                    self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
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
        pcf = self.read_pcffile()
        for i in range(0, len(self.system.list_atoms_m2)):
            for j in range(0, len(self.system.list_atoms_m2[i])):
                curr_mod = []
                for k in range(0, 3):
                    curr_mod.append(float(pcf[int(self.system.list_atoms_m2[i][j]) - 1][k]) / 0.52917721)
                curr_mod_charge = (
                    float(float(pcf[int(self.system.list_atoms_m2[i][j]) - 1][3])) - self.m2charges[i][j]
                )
                for k in range(0, len(self.system.list_atoms_qm)):
                    v1 = [
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][0] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][1] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][2] / 0.52917721,
                    ]
                    z1 = float(self.class_qm_job.qm_corrdata[k][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    for l in range(0, 3):
                        linkforces[int(self.system.list_atoms_qm[k]) - 1][l] += (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
                        linkforces[int(self.system.list_atoms_m2[i][j]) - 1][l] -= (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
                for k in range(0, len(self.system.list_coordinates_linkatoms)):
                    v1 = [
                        self.system.list_coordinates_linkatoms[k][0] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][1] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][2] / 0.52917721,
                    ]
                    z1 = float(self.class_qm_job.qm_corrdata[k + len(self.system.list_atoms_qm)][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    for l in range(0, 3):
                        linkforces[int(self.system.list_atoms_m1[k]) - 1][l] += (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
                        linkforces[int(self.system.list_atoms_m2[i][j]) - 1][l] -= (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
        m2count = 0
        linkstart = len(pcf) - 2 * len(list(_flatten(self.system.list_atoms_m2)))
        for i in range(0, len(self.system.list_atoms_m2)):
            for j in range(0, len(self.system.list_atoms_m2[i]) * 2):
                curr_mod = []
                for k in range(0, 3):
                    curr_mod.append(float(pcf[int(linkstart) + m2count][k]) / 0.52917721)
                curr_mod_charge = float(float(pcf[int(linkstart) + m2count][3]))
                m2count += 1
                for k in range(0, len(self.system.list_atoms_qm)):
                    v1 = [
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][0] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][1] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][2] / 0.52917721,
                    ]
                    z1 = float(self.class_qm_job.qm_corrdata[k][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    for l in range(0, 3):
                        linkforces[int(self.system.list_atoms_qm[k]) - 1][l] += (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
                for k in range(0, len(self.system.list_coordinates_linkatoms)):
                    v1 = [
                        self.system.list_coordinates_linkatoms[k][0] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][1] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][2] / 0.52917721,
                    ]
                    z1 = float(self.class_qm_job.qm_corrdata[k + len(self.system.list_atoms_qm)][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    for l in range(0, 3):
                        linkforces[int(self.system.list_atoms_m1[k]) - 1][l] += (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
        for i in range(0, len(self.system.list_coordinates_linkatoms)):
            v1 = [
                self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
            ]
            _flattened = list(_flatten(self.system.list_atoms_q1))
            v2 = [
                self.system.array_xyzq_current[int(_flattened[i]) - 1][0] / 0.52917721,
                self.system.array_xyzq_current[int(_flattened[i]) - 1][1] / 0.52917721,
                self.system.array_xyzq_current[int(_flattened[i]) - 1][2] / 0.52917721,
            ]
            v12 = np.array(v2) - np.array(v1)
            dist = np.linalg.norm(v12) / 0.71290813568205
            unit_vector_v12 = v12 / np.linalg.norm(v12)
            dist = dist * 0.5282272551
            forcecorr = self.databasecorrection("forces", "aminoacid_CACB", dist)
            for j in range(0, 3):
                linkforces[int(_flattened[i]) - 1][j] += -unit_vector_v12[j] * forcecorr * 0.5
                linkforces[int(self.system.list_atoms_m1[i]) - 1][j] += unit_vector_v12[j] * forcecorr * 0.5
        return linkforces
