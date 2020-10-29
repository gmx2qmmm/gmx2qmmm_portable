## To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will produce an orca.hess file from a g09 formchk -3 file.
# orca prints the internal FC hessian!

__author__ = "jangoetze"
__date__ = "$02-Jan-2018 14:45:17$"

import re
import math
import sys


def hes_xyz_fchk(inpname, outname):

    ifile = open(inpname)
    ofile = open(outname, "w")
    atom_map = {
        "1": "H",
        "2": "He",
        "3": "Li",
        "4": "Be",
        "5": "B",
        "6": "C",
        "7": "N",
        "8": "O",
        "9": "F",
        "10": "Ne",
        "11": "Na",
        "12": "Mg",
        "13": "Al",
        "14": "Si",
        "15": "P",
        "16": "S",
        "17": "Cl",
        "18": "Ar",
        "19": "K",
        "20": "Ca",
        "21": "Sc",
        "22": "Ti",
        "23": "V",
        "24": "Cr",
        "25": "Mn",
        "26": "Fe",
        "27": "Co",
        "28": "Ni",
        "29": "Cu",
        "30": "Zn",
        "31": "Ga",
        "32": "Ge",
        "33": "As",
        "34": "Se",
        "35": "Br",
        "36": "Kr",
        "37": "Rb",
        "38": "Sr",
        "39": "Y",
        "40": "Zr",
        "41": "Nb",
        "42": "Mo",
        "43": "Tc",
        "44": "Ru",
        "45": "Rh",
        "46": "Pd",
        "47": "Ag",
        "48": "Cd",
        "49": "In",
        "50": "Sn",
        "51": "Sb",
        "52": "Te",
        "53": "I",
        "54": "Xe",
        "55": "Cs",
        "56": "Ba",
        "57": "La",
        "58": "Ce",
        "59": "Pr",
        "60": "Nd",
        "61": "Pm",
        "62": "Sm",
        "63": "Eu",
        "64": "Gd",
        "65": "Tb",
        "66": "Dy",
        "67": "Ho",
        "68": "Er",
        "69": "Tm",
        "70": "Yb",
        "71": "Lu",
        "72": "Hf",
        "73": "Ta",
        "74": "W",
        "75": "Re",
        "76": "Os",
        "77": "Ir",
        "78": "Pt",
        "79": "Au",
        "80": "Hg",
        "81": "Tl",
        "82": "Pb",
        "83": "Bi",
        "84": "Po",
        "85": "At",
        "86": "Rn",
        "87": "Fr",
        "88": "Ra",
        "89": "Ac",
        "90": "Th",
        "91": "Pa",
        "92": "U",
        "93": "Np",
        "94": "Pu",
        "95": "Am",
        "96": "Cm",
        "97": "Bk",
        "98": "Cf",
        "99": "Es",
        "100": "Fm",
        "101": "Md",
        "102": "No",
        "103": "Lr",
        "104": "Rf",
        "105": "Db",
        "106": "Sg",
        "107": "Bh",
        "108": "Hs",
        "109": "Mt",
        "110": "Ds",
        "111": "Rg",
        "112": "Cn",
        "113": "Nh",
        "114": "Fl",
        "115": "Mc",
        "116": "Lv",
        "117": "Ts",
        "118": "Og",
    }
    n_a = 0
    list_i_a = []
    for line in ifile:
        match = re.search("Atomic numbers", line)
        if match:
            match = re.search("N=\s+(\d+)", line)
            n_a = int(match.group(1))
            for line in ifile:
                if len(list_i_a) >= n_a:
                    break
                splitline = line.split()
                list_i_a.extend(splitline)
            break
    ifile.close()
    ifile = open(inpname)
    list_i_c = []
    for line in ifile:
        match = re.search("Current cartesian coordinates", line)
        if match:
            match = re.search("N=\s+(\d+)", line)
            n_c = int(match.group(1))
            if n_c != (n_a * 3):
                print("Number of coordinates is not 3 times number of atoms. Exiting.")
                exit(1)
            for line in ifile:
                if len(list_i_c) >= n_c:
                    break
                splitline = line.split()
                for i in range(0, len(splitline)):
                    if (
                        math.sqrt(float(splitline[i]) * float(splitline[i]))
                        < 0.00000000000001
                    ):
                        splitline[i] = float(0.0)
                        continue
                    splitline[i] = float(splitline[i]) * 0.529177
                list_i_c.extend(splitline)
            break
    ifile.close()
    list_h_g09 = []
    n_h_g09 = 0
    ifile = open(inpname)
    for line in ifile:
        match = re.search("Cartesian Force Constants", line)
        if match:
            match = re.search("N=\s+(\d+)", line)
            n_h_g09 = int(match.group(1))
            break
    if n_h_g09 == 0:
        print("Did not find a proper Cartesian Force Constants entry. Exiting.")
        exit(1)
    n_h_g09_curr = n_h_g09
    for line in ifile:
        matchlist = re.findall("\s+([-]*\d\.\d{8}E[+,-]\d\d)", line)
        if len(matchlist) == 0:
            print("No matches for Cartesian Force Constants!")
            exit(1)
        list_h_g09 += matchlist
        n_h_g09_curr -= 5
        if n_h_g09_curr < 1:
            break
    if n_h_g09 != len(list_h_g09):
        print("Reading of g09 fchk Cartesian Force Constants failed. Exiting.")
        exit(1)
    ofile.write(
        "\n$orca_hessian_file\n\n$act_atom\n  0\n\n$act_coord\n  0\n\n$act_energy\n        0.000000\n\n$hessian\n"
    )
    n_dof = n_a * 3
    list_h_orca = []
    # initlist
    for i in range(0, n_dof):
        for j in range(0, n_dof):
            entry = 0.0
            list_h_orca.append(entry)
    counter = 0
    for i in range(0, n_dof):
        for j in range(0, i + 1):
            list_h_orca[i * n_dof + j] = list_h_g09[counter]
            list_h_orca[j * n_dof + i] = list_h_g09[counter]
            counter += 1
    ofile.write(str(n_dof) + "\n")
    rem_cols = n_dof
    curr_col = 0
    while rem_cols > 0:
        col_limit = rem_cols
        if rem_cols >= 5:
            col_limit = 5
        colstring = "  "
        for i in range(0, col_limit):
            colstring += "              {:>5}".format(curr_col + i)
        ofile.write(colstring + "\n")
        for i in range(0, n_dof):
            colstring = "{:>5}   ".format(i)
            for j in range(0, col_limit):
                colstring += "  {:>17.10E}".format(
                    float(list_h_orca[i * n_dof + curr_col + j])
                )
            ofile.write(colstring + "\n")
        curr_col += 5
        rem_cols -= 5
    ofile.write("\n$vibrational_frequencies\n" + str(n_dof) + "\n")
    ifile.close()
    ifile = open(inpname)
    n_v_entries = 0
    for line in ifile:
        match = re.search("Vib-E2", line)
        if match:
            match = re.search("N=\s+(\d+)", line)
            n_v_entries = int(match.group(1))
            break
    if n_v_entries == 0:
        print("Did not find a proper Vib-E2 entry. Exiting.")
        exit(1)
    n_v = (n_a * 3) - 6
    rem_v = n_v_entries
    list_v = []
    while rem_v > 0:
        v_limit = rem_v
        if rem_v > 5:
            v_limit = 5
        for line in ifile:
            matchlist = re.findall("\s+([-]*\d\.\d{8}E[+,-]\d\d)", line)
            for i in range(0, v_limit):
                list_v.append(float(matchlist[i]))
            break
        rem_v -= 5
    if n_v_entries != len(list_v):
        print("Reading of g09 fchk Vib-E2 failed. Exiting.")
        exit(1)
    for i in range(0, n_dof):
        colstring = "{:>5}   ".format(i)
        if i < 6:
            entry = 0.0
        else:
            entry = float(list_v[i - 6])
        colstring += "{:>13.6f}   ".format(entry)
        ofile.write(colstring + "\n")
    ofile.write("\n$normal_modes\n" + str(n_dof) + " " + str(n_dof) + "\n")
    ifile.close()
    ifile = open(inpname)
    n_nm = 0
    for line in ifile:
        match = re.search("Vib-Modes", line)
        if match:
            match = re.search("N=\s+(\d+)", line)
            n_nm = int(match.group(1))
            break
    if n_nm == 0:
        print("Did not find a proper Vib-Modes entry. Exiting.")
        exit(1)
    list_nm_g09 = []
    list_nm_orca = []
    # initlist
    for i in range(0, n_dof):
        for j in range(0, n_dof):
            entry = 0.0000000000e00
            list_nm_orca.append(entry)
    n_nm_curr = n_nm
    for line in ifile:
        matchlist = re.findall("\s+([-]*\d\.\d{8}E[+,-]\d\d)", line)
        if len(matchlist) == 0:
            print("No matches for Vib-Modes!")
            exit(1)
        list_nm_g09 += matchlist
        n_nm_curr -= 5
        if n_nm_curr < 1:
            break
    if n_nm != len(list_nm_g09):
        print("Reading of g09 fchk Vib-Modes failed. Exiting.")
        exit(1)
    for i in range(0, n_dof):
        for j in range(0, n_dof):
            if j < 6:
                list_nm_orca[i * n_dof + j] = 0.0000000000e00
            else:
                list_nm_orca[i * n_dof + j] = list_nm_g09[(n_dof * (j - 6)) + i]
    rem_cols = n_dof
    curr_col = 0
    while rem_cols > 0:
        col_limit = rem_cols
        if rem_cols >= 5:
            col_limit = 5
        colstring = "  "
        for i in range(0, col_limit):
            colstring += "              {:>5}".format(curr_col + i)
        ofile.write(colstring + "\n")
        for i in range(0, n_dof):
            colstring = "{:>5}   ".format(i)
            for j in range(0, col_limit):
                colstring += "  {:>17.10E}".format(
                    float(list_nm_orca[i * n_dof + curr_col + j])
                )
            ofile.write(colstring + "\n")
        curr_col += 5
        rem_cols -= 5
    ofile.write(
        "\n#\n# The atoms: label  mass x y z (in bohrs)\n#\n$atoms\n" + str(n_a) + "\n"
    )
    n_mass = 0
    ifile.close()
    ifile = open(inpname)
    for line in ifile:
        match = re.search("Vib-AtMass", line)
        if match:
            match = re.search("N=\s+(\d+)", line)
            n_mass = int(match.group(1))
            break
    if n_mass == 0:
        print("Did not find a proper Vib-AtMass entry. Exiting.")
        exit(1)
    n_mass_curr = n_mass
    list_mass = []
    for line in ifile:
        matchlist = re.findall("\s+([-]*\d\.\d{8}E[+,-]\d\d)", line)
        if len(matchlist) == 0:
            print("No matches for Vib-AtMass!")
            exit(1)
        list_mass += matchlist
        n_mass_curr -= 5
        if n_mass_curr < 1:
            break
    if n_mass != len(list_mass):
        print("Reading of g09 fchk Vib-AtMass failed. Exiting.")
        exit(1)
    for i in range(0, int(n_a)):
        reformatted = " {:<2} {:>11.5f}  {:>18.12f} {:>18.12f} {:>18.12f}\n".format(
            atom_map[list_i_a[i]],
            float(list_mass[i]),
            float(list_i_c[3 * i + 0]) / 0.529177,
            float(list_i_c[3 * i + 1]) / 0.529177,
            float(list_i_c[3 * i + 2]) / 0.529177,
        )
        ofile.write(reformatted)
    ofile.write(
        "\n$actual_temperature\n  0.000000\n\n$dipole_derivatives\n" + str(n_dof) + "\n"
    )
    n_dip = 0
    ifile.close()
    ifile = open(inpname)
    for line in ifile:
        match = re.search("Dipole Derivatives", line)
        if match:
            match = re.search("N=\s+(\d+)", line)
            n_dip = int(match.group(1))
            break
    if n_dip == 0:
        print("Did not find a proper Dipole Derivatives entry. Exiting.")
        exit(1)
    n_dip_curr = n_dip
    list_dip = []
    for line in ifile:
        matchlist = re.findall("\s+([-]*\d\.\d{8}E[+,-]\d\d)", line)
        if len(matchlist) == 0:
            print("No matches for Dipole Derivatives!")
            exit(1)
        list_dip += matchlist
        n_dip_curr -= 5
        if n_dip_curr < 1:
            break
    if n_dip != len(list_dip):
        print("Reading of g09 fchk Dipole Derivatives failed. Exiting.")
        exit(1)
    for i in range(0, int(n_dip) / 3):
        reformatted = "{:>20.10E}{:>20.10E}{:>20.10E}\n".format(
            float(list_dip[i * 3 + 0]),
            float(list_dip[i * 3 + 1]),
            float(list_dip[i * 3 + 2]),
        )
        ofile.write(reformatted)
    ofile.write(
        "\n#\n# The IR spectrum\n#  wavenumber T**2 TX TY  TY\n#\n$ir_spectrum\n"
        + str(n_dof)
        + "\n"
    )
    for i in range(0, 6):
        ofile.write("      0.00       0.0000       0.0000       0.0000       0.0000\n")
    for i in range(0, n_v):
        reformatted = "{:>10.2f}{:>13.4f}{:>13}{:>13}{:>13}".format(
            list_v[i], list_v[3 * n_v + i], "tbd", "tbd", "tbd"
        )
        ofile.write(reformatted + "\n")
    ofile.write("\n\n$end\n\n")
    ifile.close()
    ofile.close()


if __name__ == "__main__":
    print(hes_xyz_fchk(sys.argv[1], sys.argv[2]))
