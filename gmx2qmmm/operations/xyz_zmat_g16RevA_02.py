# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will produce a gaussian zmatrix file from an orca xyz file.

__author__ = "jangoetze"
__date__ = "$02-Jan-2018 14:45:17$"

import os
import math
import sys

import numpy as np


def make_distmat(coords):  # unclustered coords
    distmat = []
    for i in range(0, len(coords) / 3):
        curr_coords = []
        distmatline = []
        for j in range(0, 3):
            curr_coords.append(float(coords[i * 3 + j]))
        for j in range(0, len(coords) / 3):
            curr_coords2 = []
            for k in range(0, 3):
                curr_coords2.append(float(coords[j * 3 + k]))
            dist = np.array(curr_coords) - np.array(curr_coords2)
            distmatline.append(float(np.linalg.norm(dist)))
        distmat.append(distmatline)
    return distmat


def normvec(vec):
    return np.array(vec) / np.linalg.norm(vec)


def make_angle(a, b, c, coords):  # unclustered coords
    coords_a = []
    coords_b = []
    coords_c = []
    for i in range(0, 3):
        coords_a.append(float(coords[a * 3 + i]))
        coords_b.append(float(coords[b * 3 + i]))
        coords_c.append(float(coords[c * 3 + i]))
    ab = np.array(coords_a) - np.array(coords_b)
    cb = np.array(coords_c) - np.array(coords_b)
    if np.linalg.norm(ab) == 0.0 or np.linalg.norm(cb) == 0.0:
        print("Atoms with identical coordinates found. Exiting.")
        exit(1)
    angle = 180.0 * np.arccos(np.dot(normvec(ab), normvec(cb))) / math.pi
    return angle


def make_torsion(a, b, c, d, coords):  # unclustered coords
    coords_a = []
    coords_b = []
    coords_c = []
    coords_d = []
    for i in range(0, 3):
        coords_a.append(float(coords[a * 3 + i]))
        coords_b.append(float(coords[b * 3 + i]))
        coords_c.append(float(coords[c * 3 + i]))
        coords_d.append(float(coords[d * 3 + i]))
    ab = np.array(coords_a) - np.array(coords_b)
    bc = np.array(coords_b) - np.array(coords_c)
    cd = np.array(coords_c) - np.array(coords_d)
    if np.linalg.norm(ab) == 0.0 or np.linalg.norm(bc) == 0.0 or np.linalg.norm(cd) == 0.0:
        print("Atoms with identical coordinates found. Exiting.")
        exit(1)
    normab = ab / np.linalg.norm(ab)
    normbc = bc / np.linalg.norm(bc)
    normcd = cd / np.linalg.norm(cd)
    torsion = math.atan2(
        np.dot(np.cross(np.cross(normab, normbc), np.cross(normbc, normcd)), normbc),
        np.dot(np.cross(normab, normbc), np.cross(normbc, normcd)),
    )
    proper_torsion = -180.0 * torsion / math.pi
    return proper_torsion


def getbonded(atoms, current, distmat):
    mindist = float(-1.0)
    curr_min = -1
    altmindist = float(-1.0)
    altmin = -1
    proper = False
    for i in range(0, current):
        if i == current:
            continue
        if atoms[i] == "H" and (
            float(distmat[current][i]) < altmindist or altmindist < 0.0
        ):
            altmin == i
            altmindist = float(distmat[current][i])
        if atoms[i] != "H" and (float(distmat[current][i]) < mindist or mindist < 0.0):
            proper = True
            curr_min = i
            mindist = float(distmat[current][i])
    if proper:
        return curr_min
    else:
        return altmin


def getangled(atoms, coords, current, bonded, distmat):
    mindist = float(-1.0)
    curr_min = -1
    altmindist = float(-1.0)
    altmin = -1
    proper = False
    for i in range(0, current):
        if i == bonded:
            continue
        if atoms[i] == "H" and (
            float(distmat[bonded][i]) < altmindist or altmindist < 0.0
        ):
            abc = make_angle(current, bonded, i, coords)
            if abc < 0.001 or abc > 179.999:
                continue
            altmin == i
            altmindist = float(distmat[bonded][i])
        if atoms[i] != "H" and (float(distmat[bonded][i]) < mindist or mindist < 0.0):
            abc = make_angle(current, bonded, i, coords)
            if abc < 0.001 or abc > 179.999:
                continue
            proper = True
            curr_min = i
            mindist = float(distmat[bonded][i])
    if proper:
        if curr_min == -1:
            print("Only collinear atoms encountered during zmat construction. Exiting.")
            exit(1)
        else:
            return curr_min
    else:
        if altmin == -1:
            print("Only collinear atoms encountered during zmat construction. Exiting.")
            exit(1)
        return altmin


def gettorsioned(atoms, coords, current, bonded, angled, distmat):
    mindist = float(-1.0)
    curr_min = -1
    altmindist = float(-1.0)
    altmin = -1
    proper = False
    for i in range(0, current):
        if i == bonded or i == angled:
            continue
        if atoms[i] == "H" and (
            float(distmat[angled][i]) < altmindist or altmindist < 0.0
        ):
            bcd = make_angle(bonded, angled, i, coords)
            if bcd < 0.001 or bcd > 179.999:
                continue
            altmin == i
            altmindist = float(distmat[angled][i])
        if atoms[i] != "H" and (float(distmat[angled][i]) < mindist or mindist < 0.0):
            bcd = make_angle(bonded, angled, i, coords)
            if bcd < 0.001 or bcd > 179.999:
                continue
            proper = True
            curr_min = i
            mindist = float(distmat[angled][i])
    if proper:
        if curr_min == -1:
            print("Only collinear atoms encountered during zmat construction. Exiting.")
            exit(1)
        else:
            return curr_min
    else:
        if altmin == -1:
            print("Only collinear atoms encountered during zmat construction. Exiting.")
            exit(1)
        return altmin


def xyz_zmat_g16RevA02(inpfile, outfile):
    basedir = os.path.dirname(os.path.abspath(__file__))

    extrxyz = imp.load_source(
        "orcaextractors", str(basedir + "/operations/geo_from_xyz.py")
    )
    coords = extrxyz.geo_from_xyz(inpfile, "A")
    extratoms = imp.load_source(
        "orcaextractors", str(basedir + "/operations/atoms_from_xyz.py")
    )
    atoms = extratoms.atoms_from_xyz(inpfile)

    # we have all we need, now just create the zmat file. We will keep the order, but will refer only to non-H atoms.
    # create distance matrix
    distmat = make_distmat(coords)
    # create zmat matrix:
    zmat = []
    varlist = []
    count = 0
    for entry in atoms:
        count += 1
        zmatline = ""
        if count == 1:
            zmatline += " {:<2}".format(entry)
        elif count == 2:
            zmatline += " {:<2} {:>8} {:>9}".format(entry, str(1), "b2")
            varlist.append(["b2", float(distmat[1][0])])
        elif count == 3:
            if atoms[count - 1] == "H":
                zmatline += " {:<2} {:>8} {:>9} {:>8} {:>9}".format(
                    entry, str(1), "b3", str(2), "a3"
                )
                varlist.append(["b3", float(distmat[2][0])])
                varlist.append(["a3", float(make_angle(2, 0, 1, coords))])
            else:
                zmatline += " {:<2} {:>8} {:>9} {:>8} {:>9}".format(
                    entry, str(2), "b3", str(1), "a3"
                )
                varlist.append(["b3", float(distmat[2][1])])
                varlist.append(["a3", float(make_angle(2, 1, 0, coords))])
        else:
            # print zmat
            # print varlist
            # exit(1)
            # get closest bond partner
            b = getbonded(atoms, count - 1, distmat)
            bondvar = "b" + str(count)
            varlist.append([bondvar, float(distmat[count - 1][b])])
            # get closest angle partner
            a = getangled(atoms, coords, count - 1, b, distmat)
            anglevar = "a" + str(count)
            varlist.append([anglevar, float(make_angle(count - 1, b, a, coords))])
            # get closest dihedral partner
            d = gettorsioned(atoms, coords, count - 1, b, a, distmat)
            torsionvar = "d" + str(count)
            varlist.append(
                [torsionvar, float(make_torsion(count - 1, b, a, d, coords))]
            )
            zmatline = " {:<2} {:>8} {:>9} {:>8} {:>9} {:>8} {:>9}".format(
                entry, str(b + 1), bondvar, str(a + 1), anglevar, str(d + 1), torsionvar
            )
        zmat.append(zmatline)
    with open(outfile, "w") as ofile:
        for line in zmat:
            ofile.write(line + "\n")
        ofile.write("\n")
        for element in varlist:
            ofile.write(" {:<6} {:>12.6f}\n".format(element[0], float(element[1])))
        ofile.write("\n")


if __name__ == "__main__":
    xyz_zmat_g16RevA02(sys.argv[1], sys.argv[2])
