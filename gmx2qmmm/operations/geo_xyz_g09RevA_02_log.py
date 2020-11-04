# This will produce an orca.xyz file from a g09 optimisation log.

__author__ = "jangoetze"
__date__ = "$02-Jan-2018 14:45:17$"

import re
import sys


def geo_xyz_g09RevA02log(inpname, outname):
    ifile = open(inpname)
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
    linecount = 0
    droplines = True
    buffer = []
    found = False
    never_found = True
    for line in ifile:
        match = re.search("\s+Optimization completed\.", line)
        if match:
            ifile.close()
            ofile.close()
            break
        match = re.search("Input orientation", line)
        if match:
            found = True
            never_found = False
        match = re.search("Z-Matrix orientation", line)
        if match:
            found = True
            never_found = False
        match = re.search(
            "\s+\d+\s{7}\s*(\d+)\s{7}\s*\d+\s{3}\s*([-]*\d+\.\d{6})\s*([-]*\d+\.\d{6})\s*([-]*\d+\.\d{6})",
            line,
        )
        if match and found:
            if droplines:
                droplines = False
                linecount = 0
            reformatted = "  {:<2} {:>18.14f}    {:>18.14f}    {:>18.14f}".format(
                atom_map[match.group(1)],
                float(match.group(2)),
                float(match.group(3)),
                float(match.group(4)),
            )
            linecount += 1
            buffer.append("\n" + reformatted)
        elif (not match) and (not droplines):  # after a geometry
            found = False
            match = re.search(
                "\sSCF Done:\s*E\([^\s]+\)\s=\s*([-]*\d+\.\d+)\s*A\.U\.\safter", line
            )
            if match:
                ofile = open(outname, "w")
                ofile.write(
                    str(linecount)
                    + "\nCoordinates from G09 job "
                    + inpname
                    + " E "
                    + match.group(1)
                )
                ofile.write("".join(buffer))
                droplines = True
                buffer = []
                ofile.close()
    if never_found:
        ifile.close()
        ifile = open(inpname)
        for line in ifile:
            match = re.search("\s+Optimization completed\.", line)
            if match:
                ifile.close()
                ofile.close()
                break
            match = re.search("Standard orientation", line)
            if match:
                found = True
            match = re.search(
                "\s+\d+\s{7}\s*(\d+)\s{7}\s*\d+\s{3}\s*([-]*\d+\.\d{6})\s*([-]*\d+\.\d{6})\s*([-]*\d+\.\d{6})",
                line,
            )
            if match and found:
                if droplines:
                    droplines = False
                    linecount = 0
                reformatted = "  {:<2} {:>18.14f}    {:>18.14f}    {:>18.14f}".format(
                    atom_map[match.group(1)],
                    match.group(2),
                    match.group(3),
                    match.group(4),
                )
                linecount += 1
                buffer.append("\n" + reformatted)
            elif (not match) and (not droplines):  # after a geometry
                found = False
                match = re.search(
                    "\sSCF Done:\s*E\([^\s]+\)\s=\s*([-]*\d+\.\d+)\s*A\.U\.\safter",
                    line,
                )
                if match:
                    ofile = open(outname, "w")
                    ofile.write(
                        str(linecount)
                        + "\nCoordinates from G09 job "
                        + inpname
                        + " E "
                        + match.group(1)
                    )
                    ofile.write("".join(buffer))
                    droplines = True
                    buffer = []
                    ofile.close()
    if not droplines:  # in case of a full geom but no SCF Done
        ofile = open(outname, "w")
        ofile.write(
            str(linecount) + "\nCoordinates from G09 job " + inpname + " E error"
        )
        ofile.write("".join(buffer))
        ofile.close()
    ifile.close()


if __name__ == "__main__":
    geo_xyz_g09RevA02log(sys.argv[1], sys.argv[2])
