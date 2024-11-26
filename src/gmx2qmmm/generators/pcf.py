#   // INITIAL DESCRIPTION //


#   // MEATDATA //
__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$"  # During a rain storm

#   // IMPORTS //

#   Imports Of Existing Libraries
import math
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from gmx2qmmm.logging import Logger
from gmx2qmmm.generators._helper import filter_xyzq, normalized_vector, _flatten


#   // TODOS & NOTES //
#   TODO:
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class GeneratePCF():

    '''
    XX
    '''

    def __init__(self, input_dict, system, topology, directory_base) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        AJ for now I'm writing the initilization of this class as the first PCF is produced in the previous gmx2qmmm version \\
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

        self.input_dict = input_dict
        self.system = system
        self.topology = topology
        self.directory_base = directory_base

        self.pcf_filename = self.input_dict['jobname'] + ".pointcharges"


        self.make_pcf()

    def make_pcf(self):

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

        self.qm_xyzq = filter_xyzq(self.system.array_xyzq_current, self.system.list_atoms_qm, coordinates=True, charges=True)
        self.updated_chargelist = self.eliminate_and_shift_to_m1()
        new_field = self.generate_charge_shift_fieldsonly()
        self.write_new_field_to_disk_listsonly(new_field)

    def eliminate_and_shift_to_m1(self):
        '''
        ------------------------------
        EFFECT: \\
        ---------------
        removes qm charges and shifts them to m1 atoms
        ------------------------------
        INPUT: \\
        ---------------
        \\
        ------------------------------
        RETURN: \\
        ---------------
        updated_chargelist: shifts m1 charges
        ------------------------------
        '''

        updated_charges = []
        curr_charge = float(0.0)
        for element in self.qm_xyzq:
            curr_charge += float(element[3])
        count = 0
        parentcharge = float(curr_charge) - float(self.input_dict['charge'])
        for element in self.system.array_xyzq_current:
            count += 1
            if count in np.array(self.system.list_atoms_qm).astype(int):
                updated_charges.append(["QM"])
            else:
                chargeline = []
                for entry in element:
                    chargeline.append(float(entry))
                updated_charges.append(chargeline)
        if len(self.system.list_atoms_m1) != 0:
            parentcharge /= float(len(self.system.list_atoms_m1))
        count = 0
        for element in updated_charges:
            count += 1
            if count in np.array(self.system.list_atoms_m1).astype(int):
                updated_charges[count - 1][3] = float(
                    updated_charges[count - 1][3]
                ) + float(parentcharge)

        return updated_charges

    def generate_charge_shift_fieldsonly(self):

        '''
        ------------------------------
        EFFECT: \\
        ---------------
        create new pcf
        ------------------------------
        INPUT: \\
        ---------------

        ------------------------------
        RETURN: \\
        ---------------
        None
        ------------------------------
        '''

        # XX AJ no idea yet what orgfield means, but for now it contains m1 xyzq
        orgfield = filter_xyzq(self.updated_chargelist, self.system.list_atoms_m1)

        # XX AJ also no idea yet what target sum is
        target_sum = np.array([0.0, 0.0, 0.0])

        for element in self.qm_xyzq:
            target_sum += np.array(
                self.sum_pcf_tm_nofile(orgfield, element[0], element[1], element[2])
            )
        # self.updated_chargelist = pcf, target_sum same
        m1coordsq = []
        m2coordsqlist = []
        dispcharge = 0.1944  # for 0.023333333333 shifted charge on three M2 atoms; should work in any case unless distances become HUGE
        disp_charge_vec = []
        count = 0
        for element in self.system.list_atoms_m2:
            m1 = self.updated_chargelist[int(self.system.list_atoms_m1[count]) - 1].copy()
            m1coordsq.append(m1)
            m2coordsqthing = self.get_m2vec_fieldsonly(element, self.updated_chargelist)
            # XX AJ removed m2coordsq from previous version as it seems unnecessary, if it causes problems, look back into it

            if m1 == [] or m2coordsqthing == []:
                continue

            for entry in m2coordsqthing:
                entry[3] = m1[3] / float(len(m2coordsqthing))

            effect_dispcharge = float(dispcharge) * (
                m1[3] / float(len(element)) / 0.023333333333
            )

            if math.sqrt(float(effect_dispcharge) * float(effect_dispcharge)) > 0.5:
                if float(effect_dispcharge) < 0.0:
                    effect_dispcharge = -0.5
                else:
                    effect_dispcharge = 0.5
            disp_charge_vec.append(effect_dispcharge)

            count += 1
            m2coordsqlist.append(m2coordsqthing)

        disp = 0.1
        disp_vec = [disp for _ in range(len(self.system.list_atoms_m2)) for _ in range(len(self.system.list_atoms_m2[_]))]

        corr_charge_list = self.create_corr_charges(
            m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, self.system.list_atoms_m2
        )

        new_field = np.array(self.make_new_field(m2coordsqlist, corr_charge_list))
        curr_sum = np.array([0.0, 0.0, 0.0])
        for element in self.qm_xyzq:
            curr_sum += np.array(
                self.sum_pcf_tm_nofile(new_field, element[0], element[1], element[2])
            )

        curr_delta = math.sqrt(
            (curr_sum[0] - target_sum[0]) * (curr_sum[0] - target_sum[0])
            + (curr_sum[1] - target_sum[1]) * (curr_sum[1] - target_sum[1])
            + (curr_sum[2] - target_sum[2]) * (curr_sum[2] - target_sum[2])
        )
        max_change_disp = disp
        change_disp = []
        change_disp = np.array(disp_vec)
        startfac = 10.0
        v_opt_tracker = []
        for element in disp_vec:
            changestat = False
            point_vals = [0.0, 0.0, changestat]
            v_opt_tracker.append(point_vals)
        iter_count = 0
        strike = 0
        # checkpoint AJ XX temp
        while max_change_disp > 0.0000001 and iter_count < 200:

            max_change_disp = 0.0
            iter_count += 1
            count = 0
            mod_disp_vec = np.array(disp_vec)
            for i in range(0, len(self.system.list_atoms_m2)):
                for j in range(0, len(self.system.list_atoms_m2[i])):
                    if not v_opt_tracker[count][2]:
                        change_disp[count] /= startfac
                        v_opt_tracker[count][2] = True
                    if change_disp[count] <= 0.00001:
                        count += 1
                        continue
                    new_disp_vec = np.array(disp_vec)
                    new_disp_vec[count] += change_disp[count]
                    curr_corr_charge_list = self.create_corr_charges(
                        m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, self.system.list_atoms_m2
                    )
                    curr_field = np.array(
                        self.make_new_field(m2coordsqlist, curr_corr_charge_list)
                    )
                    new_sum = np.array([0.0, 0.0, 0.0])
                    for element in self.qm_xyzq:
                        new_sum += np.array(
                            self.sum_pcf_tm_nofile(
                                curr_field, element[0], element[1], element[2]
                            )
                        )
                    v_opt_tracker[count][0] = math.sqrt(
                        (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
                        + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
                        + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
                    )
                    new_disp_vec = np.array(disp_vec)
                    new_disp_vec[count] -= change_disp[count]
                    curr_corr_charge_list = self.create_corr_charges(
                        m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, self.system.list_atoms_m2
                    )
                    curr_field = np.array(
                        self.make_new_field(m2coordsqlist, curr_corr_charge_list)
                    )
                    new_sum = np.array([0.0, 0.0, 0.0])
                    for element in self.qm_xyzq:
                        new_sum += np.array(
                            self.sum_pcf_tm_nofile(
                                curr_field, element[0], element[1], element[2]
                            )
                        )
                    v_opt_tracker[count][1] = math.sqrt(
                        (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
                        + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
                        + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
                    )
                    if (float(v_opt_tracker[count][1]) < float(curr_delta)) and (
                        float(v_opt_tracker[count][1]) < float(v_opt_tracker[count][0])
                    ):
                        mod_disp_vec[count] -= change_disp[count]
                        if (
                            math.sqrt(mod_disp_vec[count] * mod_disp_vec[count])
                        ) < 0.000001:
                            mod_disp_vec[count] = 0.000001
                        v_opt_tracker[count][2] = True
                    elif (float(v_opt_tracker[count][0]) < float(curr_delta)) and (
                        float(v_opt_tracker[count][0]) < float(v_opt_tracker[count][1])
                    ):
                        mod_disp_vec[count] += change_disp[count]
                        if (
                            math.sqrt(mod_disp_vec[count] * mod_disp_vec[count])
                        ) < 0.000001:
                            mod_disp_vec[count] = 0.000001
                        v_opt_tracker[count][2] = True
                    else:
                        v_opt_tracker[count][2] = False
                    count += 1
            for element in mod_disp_vec:
                if element > max_change_disp:
                    max_change_disp = element
            new_disp_vec = np.array(mod_disp_vec)
            curr_corr_charge_list = self.create_corr_charges(
                m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, self.system.list_atoms_m2
            )
            curr_field = np.array(self.make_new_field(m2coordsqlist, curr_corr_charge_list))
            new_sum = np.array([0.0, 0.0, 0.0])
            for element in self.qm_xyzq:
                new_sum += np.array(
                    self.sum_pcf_tm_nofile(
                        curr_field, element[0], element[1], element[2]
                    )
                )
            new_delta = math.sqrt(
                (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
                + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
                + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
            )
            if new_delta < curr_delta:
                disp_vec = new_disp_vec
                curr_delta = new_delta
                curr_sum = new_sum
                strike = 0
            elif new_delta == curr_delta:
                disp_vec = new_disp_vec
                curr_delta = new_delta
                curr_sum = new_sum
                strike += 1
                if strike > 3:
                    max_change_disp = 0.0000001
                    break
            else:
                count = 0
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        change_disp[count] /= startfac
                        count += 1
            if max_change_disp < 0.0000001:
                break
        corr_charge_list = self.create_corr_charges(
            m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, self.system.list_atoms_m2
        )
        new_field = np.array(self.make_new_field(m2coordsqlist, corr_charge_list))
        new_sum = np.array([0.0, 0.0, 0.0])
        for element in self.qm_xyzq:
            new_sum += np.array(
                self.sum_pcf_tm_nofile(new_field, element[0], element[1], element[2])
            )
        new_delta = math.sqrt(
            (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
            + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
            + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
        )

        return new_field

    def sum_pcf_tm_nofile(self, inp, x, y, z):

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

        base = [float(x), float(y), float(z)]
        sumvec = np.array([0.0, 0.0, 0.0])
        for element in inp:
            if element[0] == "QM":
                continue
            new_vec = [float(element[0]), float(element[1]), float(element[2])]
            new_charge = float(element[3])
            distvec = np.array(new_vec) - np.array(base)
            qdistvec = np.array(distvec) * new_charge
            sumvec += qdistvec
        return sumvec

    def get_m2vec_fieldsonly(self, m2entry, pcf):

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

        '''extract m2 xyzq from pcf'''
        m2coordsq = []
        for i in range(0, len(m2entry)):
            m2coordsq.append(pcf[int(m2entry[i] - 1)].copy())
        return m2coordsq

    def create_corr_charges(self, m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, m2_nolist):
        '''
        ------------------------------ \\
        EFFECT: \\
        Create and normalize a list of displacement charges. \\
        --------------- \\
        INPUT: \\
        m1coordsq: float \\
        m2coordsqlist: list \\
            List of coordinates for the second masses. \\
        disp_vec: list or np.array \\
            List of displacement vectors. \\
        disp_charge_vec: list \\
            List of displacement charges. \\
        m2_nolist: list \\
            List of lists specifying the number of m2 atoms for each m1. \\
        --------------- \\
        RETURN: \\
        np.array \\
            Normalized vector of displacement charges. \\
        --------------- \\
        '''
        corr_charge_list = []
        count = 0
        for i in range(0, len(m2_nolist)):
            short_disp_vec = []
            for j in range(0, len(m2_nolist[i])):
                short_disp_vec.extend([disp_vec[count + j]])
            charge_list = self.write_disp_charges(
                m1coordsq[i], m2coordsqlist[i], short_disp_vec, disp_charge_vec[i]
            )
            corr_charge_list.extend(charge_list)
            count += len(m2_nolist[i])
        corr_charge = np.array(corr_charge_list).reshape(-1, 4)
        return corr_charge

    def write_disp_charges(self, m1, m2coordsqlist, disp_vec, dispcharge):

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

        charge_list = []
        for i in range(0, len(m2coordsqlist)):
            ab_vec = np.array(
                [
                    float(m2coordsqlist[i][0]),
                    float(m2coordsqlist[i][1]),
                    float(m2coordsqlist[i][2])
                ]
            ) - np.array([float(m1[0]), float(m1[1]), float(m1[2])])

            length_ab = np.linalg.norm(ab_vec)
            normvec = normalized_vector(ab_vec)
            ab_long = np.array(normvec) * (length_ab + disp_vec[i])
            ab_short = np.array(normvec) * (length_ab - disp_vec[i])
            new_dispcharge = 0.0
            if disp_vec[i] >= 0.000001:
                new_dispcharge = 0.5 - (0.5 / (2.0 * disp_vec[i] + 1.0))
            if dispcharge > 0.0:
                chargeminus = [
                    float(m1[0]) + float(ab_long[0]),
                    float(m1[1]) + float(ab_long[1]),
                    float(m1[2]) + float(ab_long[2]),
                    float(-1.0 * new_dispcharge),
                ]
                chargeplus = [
                    float(m1[0]) + float(ab_short[0]),
                    float(m1[1]) + float(ab_short[1]),
                    float(m1[2]) + float(ab_short[2]),
                    float(new_dispcharge),
                ]
            else:
                chargeplus = [
                    float(m1[0]) + float(ab_long[0]),
                    float(m1[1]) + float(ab_long[1]),
                    float(m1[2]) + float(ab_long[2]),
                    float(new_dispcharge),
                ]
                chargeminus = [
                    float(m1[0]) + float(ab_short[0]),
                    float(m1[1]) + float(ab_short[1]),
                    float(m1[2]) + float(ab_short[2]),
                    float(-1.0 * new_dispcharge),
                ]
            charge_list.append(chargeminus)
            charge_list.append(chargeplus)
        return charge_list

    def make_new_field(self, m2coordsqlist,corr_charge_list):

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

        '''combine m2 and corr charges to a new list'''
        new_field=[]
        for element in m2coordsqlist:
            for entry in element:
                new_field.append(entry)
        for element in corr_charge_list:
            new_field.append(element)
        return new_field

    def write_new_field_to_disk_listsonly(self, new_field):
        '''
        ------------------------------
        EFFECT: \\
        ---------------
        write point charge field file
        ------------------------------
        INPUT: \\
        ---------------
        inp: list, list, xyzq for all atoms, 'QM' for qm atoms
        ofilename: string, name for pcf file
        new_field: 2d array, xyzq of new charges
        getlist: list, m1 atoms
        m2_nolist: list, m2 atoms
        ------------------------------
        RETURN: \\
        ---------------
        None
        ------------------------------
        '''

        with open(self.pcf_filename, "w") as ofile:
            count = 0
            list_atoms_m2 = np.array(self.system.list_atoms_m2).reshape(-1)
            m2count = 0
            for element in self.updated_chargelist:
                count += 1
                if int(count) in np.array(list(_flatten(list_atoms_m2))).astype(int):
                    if len(element) != 4:
                        print("Line " + str(count) + " does not contain data. Exiting.")
                        exit(1)
                    ofile.write(
                        "{:<.10f} {:<.10f} {:<.10f} {:<.10f}\n".format(
                            float(element[0]),
                            float(element[1]),
                            float(element[2]),
                            float(new_field[m2count][3]) + float(element[3]),
                        )
                    )
                    m2count += 1
                    continue
                if int(count) in np.array(self.system.list_atoms_m1).astype(int):
                    ofile.write("QM\n")
                    continue
                else:
                    for i in range(0, len(element)):
                        ofile.write(str(element[i]))
                        if i != len(element) - 1:
                            ofile.write(" ")
                    ofile.write("\n")
            for i in range(len(list(_flatten(list_atoms_m2))), len(new_field)):
                ofile.write(
                    "{:<.10f} {:<.10f} {:<.10f} {:<.10f}\n".format(
                        float(new_field[i][0]),
                        float(new_field[i][1]),
                        float(new_field[i][2]),
                        float(new_field[i][3]),
                    )
                )
            ofile.write("$end\n")
