# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will read an orca.hess file and convert atom types to masses, reporting them.

# for some elements, the most stable isotope was chosen (all from Nuclear Physics A 729: 3-128): 98Tc, 147Pm, 209Po, 210At, 222Rn, 223Fr, 226Ra, 227Ac, 237Np, 244Pu, 243Am, 247Cm, 247Bk, 251Cf, 252Es, 257Fm, 258Md, 259No, 266Lr, 267Rf, 268Db, 269Sg, 270Bh, 277Hs, 278Mt, 281Ds, 282Rg, 285Cn, 286Nh, 289Fl, 289Mc, 293Lv, 294Ts, 295Og
# source: IUPAC PTE as of 29-11-2016, IUPAC website

__author__ = "jangoetze"
__date__ = "$02-Jan-2018 14:45:17$"

import re
import sys


def mass_from_hes(inpname):
    ifile = open(inpname)
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
    n_a = 0
    list_m = []
    for line in ifile:
        match = re.search("\$atoms", line)
        if match:
            for line in ifile:
                n_a = int(line)
                break
            for line in ifile:
                match = re.search(
                    "\s+(\w+)\s+\d+\.\d{5}\s+[-]*\d+\.\d{12}\s+[-]*\d+\.\d{12}\s+[-]*\d+\.\d{12}",
                    line,
                )
                list_m.append(mass_map[match.group(1)])
                if len(list_m) == n_a:
                    break
            break
    ifile.close()
    if len(list_m) != n_a:
        print("Did not find the same number of atoms as indicated by number of atoms at beginning of geometry!")
        exit(1)
    return list_m


if __name__ == "__main__":
    mass_from_hes(sys.argv[1])
