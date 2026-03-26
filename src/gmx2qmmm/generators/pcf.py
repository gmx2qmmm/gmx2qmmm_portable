"""Generate point charge fields (PCF) for QM/MM calculations"""

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$"  # During a rain storm

import math
from collections.abc import Mapping
from functools import cached_property
from io import TextIOBase
from typing import List, Self

import numpy as np

from gmx2qmmm.generators._helper import filter_xyzq, normalized_vector, _flatten
from gmx2qmmm.generators.system import SystemInfo
from gmx2qmmm.generators.types import Point, Points, PointChargeField, Vector, Vectors3D


class GeneratePCF:
    """ "Deprecated class for PCF generation

    Use :class:`PCF` instead.
    """

    def __init__(self, input_dict, system, topology, work_dir) -> None:
        """
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
        """

        self.input_dict = input_dict
        self.system = system
        self.topology = topology
        self.work_dir = work_dir

        self.pcf_filename = self.work_dir / (
            self.input_dict["jobname"] + ".pointcharges"
        )

        self.make_pcf()

    def make_pcf(self):
        """
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
        """

        self.qm_xyzq = filter_xyzq(
            self.system.array_xyzq_current,
            self.system.list_atoms_qm,
            coordinates=True,
            charges=True,
        )
        self.updated_chargelist = self.eliminate_and_shift_to_m1()
        new_field = self.generate_charge_shift_fieldsonly()
        self.write_new_field_to_disk_listsonly(new_field)

    def eliminate_and_shift_to_m1(self):
        """
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
        """

        updated_charges = []
        curr_charge = float(0.0)
        for element in self.qm_xyzq:
            curr_charge += float(element[3])
        count = 0
        parentcharge = float(curr_charge) - float(self.input_dict["charge"])
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
        """
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
        """

        # XX AJ no idea yet what orgfield means, but for now it contains m1 xyzq
        orgfield = filter_xyzq(self.updated_chargelist, self.system.list_atoms_m1)

        # XX AJ also no idea yet what target sum is
        target_sum = np.array([0.0, 0.0, 0.0])

        for element in self.qm_xyzq:
            target_sum += np.array(
                self.sum_pcf_tm_nofile(orgfield, element[0], element[1], element[2])
            )

        # For later inspection
        self.target_sum = target_sum

        # self.updated_chargelist = pcf, target_sum same
        m1coordsq = []
        m2coordsqlist = []
        dispcharge = 0.1944  # for 0.023333333333 shifted charge on three M2 atoms; should work in any case unless distances become HUGE
        disp_charge_vec = []
        count = 0
        for element in self.system.list_atoms_m2:
            m1 = self.updated_chargelist[
                int(self.system.list_atoms_m1[count]) - 1
            ].copy()
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
        disp_vec = [
            disp
            for _ in range(len(self.system.list_atoms_m2))
            for _ in range(len(self.system.list_atoms_m2[_]))
        ]

        corr_charge_list = self.create_corr_charges(
            m1coordsq,
            m2coordsqlist,
            disp_vec,
            disp_charge_vec,
            self.system.list_atoms_m2,
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
                        m1coordsq,
                        m2coordsqlist,
                        new_disp_vec,
                        disp_charge_vec,
                        self.system.list_atoms_m2,
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
                        m1coordsq,
                        m2coordsqlist,
                        new_disp_vec,
                        disp_charge_vec,
                        self.system.list_atoms_m2,
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
                m1coordsq,
                m2coordsqlist,
                new_disp_vec,
                disp_charge_vec,
                self.system.list_atoms_m2,
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

        self.iterations_ = iter_count

        corr_charge_list = self.create_corr_charges(
            m1coordsq,
            m2coordsqlist,
            disp_vec,
            disp_charge_vec,
            self.system.list_atoms_m2,
        )
        new_field = np.array(self.make_new_field(m2coordsqlist, corr_charge_list))
        new_sum = np.array([0.0, 0.0, 0.0])
        self.new_sum = new_sum

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
        """
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
        """

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
        """
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
        """

        """extract m2 xyzq from pcf"""
        m2coordsq = []
        for i in range(0, len(m2entry)):
            m2coordsq.append(pcf[int(m2entry[i] - 1)].copy())
        return m2coordsq

    def create_corr_charges(
        self, m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, m2_nolist
    ):
        """
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
        """
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
        """
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
        """

        charge_list = []
        for i in range(0, len(m2coordsqlist)):
            ab_vec = np.array(
                [
                    float(m2coordsqlist[i][0]),
                    float(m2coordsqlist[i][1]),
                    float(m2coordsqlist[i][2]),
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

    def make_new_field(self, m2coordsqlist, corr_charge_list):
        """
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
        """

        """combine m2 and corr charges to a new list"""
        new_field = []
        for element in m2coordsqlist:
            for entry in element:
                new_field.append(entry)
        for element in corr_charge_list:
            new_field.append(element)
        return new_field

    def write_new_field_to_disk_listsonly(self, new_field):
        """
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
        """

        ofile = open(self.pcf_filename, "w")
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
        ofile.close()


class PCFGeneratorShift:
    """Interface for point charge field (PCF) generation

    Starting from an original PCF (particle coordinates plus partial
    charges) of a full system (in the sense of including QM, M1, and M2
    atoms),
    the goal is to generate a new PCF
    for the QM/MM system under charge neutrality, where the charges of
    the QM atoms are removed and eventually redistributed to the M2
    atoms. Corrective dipol charges are added and optimised for all
    M1-M2 pairs to minimize the difference between the original and new
    PCF felt by the QM atoms.

    To this effect, the current default implementation follows the
    following procedure:
        * Remove the charges of the QM atoms and shift them to the M1
        atoms.
        * Remove the charges of the M1 atoms and shift them to the M2 atoms.
        * Create displacement charges for each M1-M2 pair,
        where the displacement charges are placed along the M1-M2 vector
        * Optimize the positions of the displacement charges to minimize
        the difference between the original and new PCF felt by the QM
        atoms
    """

    _default_parameters = {
        "initial_displacement_charge": 0.1944,
        # Corresponds to 0.0233333333333 shifted charge on three M2
        # atoms; initial guess to start optimisation
        "initial_displacement": 0.1,
        "initial_displacement_delta": 0.01,
        "displacement_modifier": 0.1,
        "reference_charge_per_m2": 0.0233333333333,
        "max_displacement_charge": 0.5,
        # Effective cap for displacement charges
        "epsilon": 1e-7,
        # Convergence criterion for optimization of displacement
        # charges; the optimization is considered converged if the maximum
        # change in displacement is less than this value
        "min_displacement": 1e-6,
        # Minimum allowed displacement after optimization of displacement charges
        "max_iterations": 200,
        # Displacement change modifier for optimization of displacement charges;
        # the change in displacement is divided by this value after each iteration
        # that failed to decrease the error to improve convergence
    }

    def __init__(
        self,
        input_field: PointChargeField,
        qm_atoms: List[int],
        m1_atoms: List[int],
        m2_atoms: List[List[int]],
        charge: float = 0,
        **parameters: float | int
    ) -> None:
        """Ïnitialize the PCF generator

        Args:
            input_field: Original point charge field of the full system,
                including QM, M1, and M2 atoms.
            qm_atoms: 0-based indices of the QM atoms in the input field.
            m1_atoms: 0-based indices of the M1 atoms in the input field.
            m2_atoms: 0-based indices of the M2 atoms for each M1 atom
                in the input field.
            charge: Total charge of the system.
            parameters: Optional parameters for the generator, which can be used
                to override the default parameters
                (see :attr:`PCFGeneratorShift._default_parameters`).
        """

        self._reset = False
        self.input_field = input_field
        self.qm_atoms = qm_atoms
        self.m1_atoms = m1_atoms
        self.m2_atoms = m2_atoms
        self.charge = charge

        self._parameters = {**self._default_parameters, **parameters}

    @classmethod
    def from_system(cls, system_info: SystemInfo, charge: float = 0, **parameters) -> Self:
        """Create a PCF instance from a SystemInfo object

        As atom indices in SystemInfo are 1-based,
        they are converted to 0-based indices for the PCF generator.

        Args:
            system_info: :class:`SystemInfo` object
            charge: Total charge of the system.
            parameters: Optional parameters for the generator, which can be used
                to override the default parameters
                (see :attr:`PCFGeneratorShift._default_parameters`).
        """

        qm_atoms = [index - 1 for index in system_info.list_atoms_qm]
        m1_atoms = [index - 1 for index in system_info.list_atoms_m1]
        m2_atoms = [[index - 1 for index in group] for group in system_info.list_atoms_m2]

        return cls(
            input_field=system_info.array_xyzq_current,
            qm_atoms=qm_atoms,
            m1_atoms=m1_atoms,
            m2_atoms=m2_atoms,
            charge=charge,
            **parameters,
        )

    @property
    def input_field(self) -> PointChargeField:
        return self._input_field

    @input_field.setter
    def input_field(self, value: PointChargeField) -> None:
        """Set the input point charge field and trigger reset of the generator"""
        if not isinstance(value, np.ndarray):
            raise TypeError("PCF must be a numpy array")
        if value.ndim != 2 or value.shape[1] != 4:
            raise ValueError("PCF must be a 2D array with shape (n_atoms, 4)")
        self._input_field = value
        self.reset()

    @property
    def qm_atoms(self) -> List[int]:
        return self._qm_atoms

    @qm_atoms.setter
    def qm_atoms(self, value: List[int]) -> None:
        """Set the QM atoms and trigger reset of the generator"""
        if not isinstance(value, list) or not all(isinstance(i, int) for i in value):
            raise TypeError("QM atoms must be a list of integers")
        self._qm_atoms = value
        self.reset()

    @property
    def m1_atoms(self) -> List[int]:
        return self._m1_atoms

    @m1_atoms.setter
    def m1_atoms(self, value: List[int]) -> None:
        """Set the M1 atoms and trigger reset of the generator"""
        if not isinstance(value, list) or not all(isinstance(i, int) for i in value):
            raise TypeError("M1 atoms must be a list of integers")
        self._m1_atoms = value
        self.reset()

    @property
    def m2_atoms(self) -> List[List[int]]:
        return self._m2_atoms

    @m2_atoms.setter
    def m2_atoms(self, value: List[List[int]]) -> None:
        """Set the M2 atoms and trigger reset of the generator"""

        if not isinstance(value, list) or not all(
            isinstance(i, list) and all(isinstance(j, int) for j in i) for i in value
        ):
            raise TypeError("M2 atoms must be a list of lists of integers")
        self._m2_atoms = value
        self.reset()

    @property
    def charge(self) -> float:
        return self._charge

    @charge.setter
    def charge(self, value: float) -> None:
        """Set the total charge of the system and trigger reset of the generator"""
        self._charge = float(value)
        self.reset()

    @property
    def qm_charges(self) -> Vector:
        return self.input_field[self.qm_atoms, 3]

    @property
    def qm_positions(self) -> Points:
        return self.input_field[self.qm_atoms, :3]

    @property
    def qm_field(self) -> PointChargeField:
        return self.input_field[self.qm_atoms]

    @property
    def m1_charges(self) -> Vector:
        return self.input_field[self.m1_atoms, 3]

    @property
    def m1_positions(self) -> Points:
        return self.input_field[self.m1_atoms, :3]

    @property
    def m1_field(self) -> PointChargeField:
        return self.input_field[self.m1_atoms]

    @cached_property
    def m2_indices(self) -> List[int]:
        return [index for group in self.m2_atoms for index in group]

    @cached_property
    def m1m2_pairs(self) -> List[int]:
        return [(m1_index, m2_index) for m1_index, m2_indices in zip(self.m1_atoms, self.m2_atoms) for m2_index in m2_indices]

    @cached_property
    def correction_indices(self) -> List[int]:
        n = self.input_field.shape[0]
        return list(range(n, n + 2 * len(self.m2_indices)))

    @cached_property
    def m2_and_correction_indices(self) -> List[int]:
        return self.m2_indices + self.correction_indices

    @property
    def m2_charges(self) -> Vector:
        return self.input_field[self.m2_indices, 3]

    @property
    def m2_positions(self) -> PointChargeField:
        return self.input_field[self.m2_indices, :3]

    @property
    def m2_field(self) -> PointChargeField:
        return self.input_field[self.m2_indices]

    @cached_property
    def annotations(self) -> List[str]:

        mapping = {}
        for index in self.qm_atoms:
            mapping[index] = "QM"
        for index in self.m1_atoms:
            mapping[index] = "M1"
        for index in self.m2_indices:
            mapping[index] = "M2"
        for index in self.correction_indices:
            mapping[index] = "CORRECTION"

        return [mapping.get(i, "OTHER") for i in range(self.input_field.shape[0] + len(self.correction_indices))]

    def reset(self) -> None:
        """Reset the generator to the original state

        No-op if the generator is already in the original state,
        otherwise clear all generated data (attributes suffixed with
        "_")
        and cached properties to ensure that the next call to
        :func:`generate` will produce the correct output based on
        possibly updated input data
        """

        if self._reset:
            return

        self.output_field_ = None
        self.qm_total_charge_ = None
        self.target_field_vector_ = None
        self.current_field_vector_ = None

        # Clear cached properties
        if hasattr(self, "m2_indices"):
            del self.m2_indices
        if hasattr(self, "m1m2_pairs"):
            del self.m1m2_pairs
        if hasattr(self, "correction_indices"):
            del self.correction_indices
        if hasattr(self, "m2_and_correction_indices"):
            del self.m2_and_correction_indices

        self._reset = True

    def generate(self) -> PointChargeField:

        self.output_field_ = np.vstack([
            self.input_field.copy(),
            np.zeros((len(self.m2_indices) * 2, 4))  # Placeholder for M2 correction charges
        ])

        self._shift_qm_charges_to_m1()
        self._shift_m1_charges_to_m2()

        return self.output_field_

    def _shift_qm_charges_to_m1(self) -> None:
        """Remove QM charges and shift them to M1 atoms"""

        self.qm_total_charge_ = np.sum(self.qm_charges)
        self.output_field_[self.qm_atoms, 3] = 0

        if len(self.m1_atoms) > 0:
            qm_parent_charge = self.qm_total_charge_ - self.charge
            charge_per_m1 = qm_parent_charge / len(self.m1_atoms)
            self.output_field_[self.m1_atoms, 3] += charge_per_m1

    def _shift_m1_charges_to_m2(self) -> None:
        """Shift M1 charges to M2 atoms using displacement charges

        Assume that the M1 charges have already been updated by
        :func:`_shift_qm_charges_to_m1`.
        """

        disp_charge = self._parameters["initial_displacement_charge"]
        ref_charge = self._parameters["reference_charge_per_m2"]
        disp = self._parameters["initial_displacement"]
        disp_delta = self._parameters["initial_displacement_delta"]
        max_disp_charge = self._parameters["max_displacement_charge"]
        eps = self._parameters["epsilon"]
        min_disp = self._parameters["min_displacement"]
        max_iter = self._parameters["max_iterations"]
        disp_modifier = self._parameters["displacement_modifier"]

        field = self.output_field_
        m1_field = field[self.m1_atoms]

        # Compute the total effective charge-weighted displacement vector
        # from the M1 field acting on the QM positions
        self.target_field_vector_ = self._effective_charge_weighted_displacement_vector(
            self.qm_positions, m1_field
        )

        # Distribute M1 charges to respective M2 atoms
        for m1, m2_indices in zip(self.m1_atoms, self.m2_atoms):
            m1_charge = field[m1, 3]
            charge_per_m2 = m1_charge / len(m2_indices)
            field[m1, 3] = 0
            # Permanently set M1 charges to 0
            field[m2_indices, 3] = charge_per_m2
            # Temporarily set M2 charges to the shifted M1 charge for
            # The optimization of displacement charges
            # Original M2 charges will be added back after optimization of displacement charges

        # Prepare displacement vectors for optimization (1 element per M1-M2 pair)
        disp_charge_vector = disp_charge * field[self.m2_indices, 3] / ref_charge
        np.clip(disp_charge_vector, -max_disp_charge, max_disp_charge, out=disp_charge_vector)
        disp_vector = np.full_like(disp_charge_vector, disp)

        # Compute two initial correction charges for each M1-M2 pair and add them to the output field
        self._update_all_correction_charges(disp_vector, disp_charge_vector)
        delta = self._objective()

        # Opt loop
        iteration = 0
        strike = 0
        max_disp_change = np.inf
        delta_disp_vector = np.full_like(disp_charge_vector, disp_delta)
        while iteration < max_iter and max_disp_change > eps:

            new_disp_vector = np.copy(disp_vector)
            new_delta_disp_vector = np.copy(delta_disp_vector)

            iteration += 1
            max_disp_change = 0

            for i, (m1, m2) in enumerate(self.m1m2_pairs):

                short_index = self.correction_indices[2 * i]
                long_index = short_index + 1
                current_short = np.copy(field[short_index])
                current_long = np.copy(field[long_index])

                # Trial displacements in both directions for the current M1-M2 pair
                self._get_correction_charges(
                    m1=field[m1, :3],
                    m2=field[m2, :3],
                    displacement=new_disp_vector[i] + new_delta_disp_vector[i],
                    charge=disp_charge_vector[i],
                    out_short=field[short_index],
                    out_long=field[long_index],
                )

                new_delta_plus = self._objective()

                self._get_correction_charges(
                    m1=field[m1, :3],
                    m2=field[m2, :3],
                    displacement=new_disp_vector[i] - new_delta_disp_vector[i],
                    charge=disp_charge_vector[i],
                    out_short=field[short_index],
                    out_long=field[long_index],
                )

                new_delta_minus = self._objective()

                if new_delta_plus < delta and new_delta_plus < new_delta_minus:
                    new_disp_vector[i] += new_delta_disp_vector[i]
                elif new_delta_minus < delta and new_delta_minus < new_delta_plus:
                    new_disp_vector[i] -= new_delta_disp_vector[i]
                else:
                    # No improvement; reduce the displacement change for the next iteration
                    new_delta_disp_vector[i] = max(new_delta_disp_vector[i] * disp_modifier, min_disp)

                field[short_index] = current_short
                field[long_index] = current_long

            # Evaluate the new configuration after iterating through all
            # M1-M2 pairs individually
            self._update_all_correction_charges(new_disp_vector, disp_charge_vector)
            new_delta = self._objective()

            if new_delta < delta:
                delta = new_delta
                self.current_field_vector_ = self.current_field_vector_trial
                max_disp_change = np.max(np.abs(new_disp_vector - disp_vector))
                disp_vector = new_disp_vector
                delta_disp_vector = new_delta_disp_vector
                strike = 0
            elif new_delta == delta:
                disp_vector = new_disp_vector
                delta_disp_vector = new_delta_disp_vector
                strike += 1
                if strike > 3:
                    break
            else:
                max_disp_change = np.inf
                delta_disp_vector *= disp_modifier
                np.clip(disp_vector, min_disp, None, out=disp_vector)

        self._update_all_correction_charges(disp_vector, disp_charge_vector)
        field[self.m2_indices, 3] += self.m2_charges
        # Adding back the original M2 charges

        self.iterations_ = iteration

        return field

    def _objective(self) -> float:
        """Calculate the objective function value for the current field configuration"""

        self.current_field_vector_trial = self._effective_charge_weighted_displacement_vector(
            self.qm_positions,
            field=self.output_field_[self.m2_and_correction_indices],
        )
        return np.linalg.norm(self.current_field_vector_trial - self.target_field_vector_)

    @staticmethod
    def _charge_weighted_displacement_vectors(
        point: Point, field: PointChargeField
    ) -> Vectors3D:
        """Calculate charge-weighted displacement vectors

        For a given point (particle coordinates) and a point charge field
        (array of particle coordinates and charges),
        calculate displacement vectors from that
        point to each charge in the field, multiplied by that charge.

        Args:
            point: Target point for which to calculate displacements
            field: Point charge field to calculate displacements from

        Returns:
            Array of shape (n_atoms, 3) containing charge-weighted
            displacement vectors from the input point to each charge
            in the field
        """

        displacements = field[:, :3] - point
        charges = field[:, 3]
        return displacements * charges[:, np.newaxis]

    def _effective_charge_weighted_displacement_vector(
        self, points: Union[Point, Points], field: PointChargeField
    ) -> Vector:
        """Calculate an effective charge-weighted displacement vector
        from multiple contributions

        For a given (set of) point(s) (particle coordinates)
        and a point charge field
        (array of particle coordinates and charges),
        calculate an effective displacement vector from
        displacement vectors of
        individual points to all charges in the field,
        multiplied by those respective charges.

        Args:
            points: Target point(s) for which to calculate displacements
            field: Point charge field to calculate displacements from

        Returns:
            Array of shape (3,) containing the total charge-weighted
            displacement vector from the input point to all charges
            in the field
        """

        def _1d_wrapper(point, field):
            return self._charge_weighted_displacement_vectors(point, field).sum(axis=0)

        if points.ndim == 1:
            return _1d_wrapper(points, field)

        return np.apply_along_axis(
            _1d_wrapper,
            1,
            points,
            field=field,
        ).sum(axis=0)

    def _update_all_correction_charges(self, disp_vector: Vector, disp_charge_vector: Vector) -> None:
        """Calculate and update the correction charges for all M1-M2 pairs

        Args:
            disp_vector: Array of shape (n_pairs,) containing the displacements for each M1-M2 pair
            disp_charge_vector: Array of shape (n_pairs,) containing the
                charge magnitudes for the correction charges of each M1-M2 pair
        """

        field = self.output_field_
        for i, (m1, m2) in enumerate(self.m1m2_pairs):
            short_index = self.correction_indices[2 * i]
            long_index = short_index + 1

            self._get_correction_charges(
                m1=field[m1, :3],
                m2=field[m2, :3],
                displacement=disp_vector[i],
                charge=disp_charge_vector[i],
                out_short=field[short_index],
                out_long=field[long_index],
            )

    def _get_correction_charges(
        self,
        m1: Point,
        m2: Point,
        displacement: float,
        charge: float,
        out_short: Optional[PointCharge] = None,
        out_long: Optional[PointCharge] = None,
    ) -> Tuple[PointCharge, PointCharge]:
        """Calculate position and magnitude of two correction charges
        for a given M1-M2 pair

        Args:
            m1: Coordinates of the M1 atom
            m2: Coordinates of the M2 atom
            displacement: Displacement along the M1-M2 vector for placing the correction charges
            charge: Charge magnitude for the correction charges.
                Positive value corresponds to a positive charge on the
                short side and negative charge on the long side, and vice versa.
            out_short: Optional output array to store the short
                correction point charge (placed at d(M1, M2) - displacement from M1)
            out_long: Optional output array to store the long
                correction point charge (placed at d(M1, M2) + displacement from M1)

        Returns:
            short and long correction point charges as 1D arrays of
            shape (4,) with elements (x, y, z, q)
        """

        if out_short is None:
            out_short = np.zeros(4)

        if out_long is None:
            out_long = np.zeros(4)

        ab = m2 - m1
        length = np.linalg.norm(ab)
        ab_normed = ab / length
        # Have to assume length > 0 (M1 and M2 atoms cannot be at the same position)

        out_short[:3] = m1 + ab_normed * (length - displacement)
        out_long[:3] = m1 + ab_normed * (length + displacement)

        new_charge = 0.5 - (0.5 / (2.0 * displacement + 1.0))

        if charge > 0.0:
            out_short[3] = new_charge
            out_long[3] = -new_charge
        else:
            out_short[3] = -new_charge
            out_long[3] = new_charge

        return out_short, out_long


def dump_field(
    fp: TextIOBase,
    field: PointChargeField,
    precision: int = 10,
    annotations: Optional[List[str]] = None,
    ) -> None:
    """Dump a point charge field to a file-like object in the format of a PCF file"""

    if annotations is not None:
        if len(annotations) != field.shape[0]:
            raise ValueError(f"Length of annotations {len(annotations)} must match the number of charges in the field {field.shape[0]}")

    charge_str = f"{{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}}"
    for i, charge in enumerate(field):
        annotation_str = ""
        if annotations is not None:
            annotation_str = f"  # {annotations[i]}"
        fp.write(charge_str.format(*charge) + annotation_str + "\n")


def load_field(fp: TextIOBase) -> PointChargeField:
    """Load a point charge field from a file-like object in the format of a PCF file"""

    charges = []
    for line in fp:
        line = line.strip()
        if not line or line.startswith("$end"):
            continue

        if line.startswith("QM"):
            # Legacy handling for "QM" lines
            line = "0.0 0.0 0.0 0.0"

        # Strip off annotation if present (anything after a "#" character)
        try:
            data, annotation = line.split("#", 1)
            annotation = annotation.strip()
        except ValueError:
            data = line
            annotation = None

        parts = data.split()
        if len(parts) != 4:
            raise ValueError(f"Invalid line in PCF file: {line}")
        try:
            parsed = [float(part) for part in parts[:4]]
        except ValueError as e:
            raise ValueError(f"Invalid numeric value in line: {line}") from e

        charges.append(parsed)

    return np.array(charges)
