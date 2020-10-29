# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will diagonalize a matrix. should be made general, currently taking the hessian from hessian_from_hes and mass weight

__author__ = "jangoetze"
__date__ = "$06-Feb-2018 12:45:17$"

import sys
import imp
import math
from numpy import linalg as LA
from numpy import array as arr
from numpy import dot as dot
from numpy import transpose as transpose
from numpy import sqrt as sqr


def make_stretch(curr_index, shaped_coords, basedir, n_a):
    import imp

    rot = imp.load_source("operations", str(basedir + "/operations/expansion_check.py"))
    stretchline = []
    for j in range(0, curr_index - 1):
        stretchline.append([0.0, 0.0, 0.0])
    unit_fwd = rot.uvec(
        arr(shaped_coords[curr_index]) - arr(shaped_coords[curr_index - 1])
    )
    unit_bwd = arr(unit_fwd) * -1.0
    ufwd = []
    ubwd = []
    for element in unit_fwd:
        ufwd.append(float(element))
    for element in unit_bwd:
        ubwd.append(float(element))
    stretchline.append(ufwd)
    stretchline.append(ubwd)
    for j in range(curr_index + 1, n_a):
        stretchline.append([0.0, 0.0, 0.0])
    return stretchline


def make_angle(curr_index, shaped_coords, basedir, n_a):
    import imp

    rot = imp.load_source("operations", str(basedir + "/operations/expansion_check.py"))
    angleline = []
    for j in range(0, curr_index - 2):
        angleline.append([0.0, 0.0, 0.0])
    unit_v1 = rot.uvec(
        arr(shaped_coords[curr_index]) - arr(shaped_coords[curr_index - 1])
    )
    unit_v2 = rot.uvec(
        arr(shaped_coords[curr_index - 2]) - arr(shaped_coords[curr_index - 1])
    )
    phi = rot.angle3d(unit_v1, unit_v2)
    angleelement = (math.cos(phi) * unit_v1 - unit_v2) / math.sin(phi)
    angleline.append([angleelement[0], angleelement[1], angleelement[2]])
    distv1 = arr(shaped_coords[curr_index]) - arr(shaped_coords[curr_index - 1])
    distv2 = arr(shaped_coords[curr_index - 2]) - arr(shaped_coords[curr_index - 1])
    sqrtvec = [
        sqr(abs(distv1[0] * distv2[0])),
        sqr(abs(distv1[1] * distv2[1])),
        sqr(abs(distv1[2] * distv2[2])),
    ]
    if distv1[0] * distv2[0] < 0.0:
        sqrtvec[0] *= -1.0
    if distv1[1] * distv2[1] < 0.0:
        sqrtvec[1] *= -1.0
    if distv1[2] * distv2[2] < 0.0:
        sqrtvec[2] *= -1.0
    angleelement = (distv1 - math.cos(phi) * distv2) * unit_v1 + unit_v2 * (
        distv2 - distv1 * math.cos(phi)
    ) / (arr(sqrtvec) * math.sin(phi))
    angleline.append([angleelement[0], angleelement[1], angleelement[2]])
    angleelement = (math.cos(phi) * unit_v2 - unit_v1) / math.sin(phi)
    angleline.append([angleelement[0], angleelement[1], angleelement[2]])
    for j in range(curr_index + 1, n_a):
        angleline.append([0.0, 0.0, 0.0])
    return angleline


def construct_bmat(coords, basedir):
    from compiler.ast import flatten

    n_a = len(coords) / 3
    shaped_coords = []
    for i in range(0, n_a):
        coordline = []
        for j in range(0, 3):
            coordline.append(coords[i * 3 + j])
        shaped_coords.append(arr(coordline))
    b_mat = []
    for i in range(1, n_a):
        if n_a < 3:
            print "Molecule too small. Exiting."
            exit(1)
        b_mat_line_stretch = make_stretch(i, shaped_coords, basedir, n_a)
        b_mat.append(flatten(b_mat_line_stretch))
        if i == 1:
            continue
        b_mat_line_angle = make_angle(i, shaped_coords, basedir, n_a)
        b_mat.append(flatten(b_mat_line_angle))
    return b_mat


def proj(u, v):
    absu = sqr(dot(u, u))
    if absu < 0.000000000000001:
        p = float(0.0) * arr(u)
    else:
        p = (dot(u, v) / dot(u, u)) * arr(u)
    return p


def newvec(mat, dim, rdm):
    rdmvec = []
    for j in range(0, int(dim)):
        rdmvec.append(float(10.0))
        if j == int(rdm):
            rdmvec[j] += float(20.0)
    scale = dot(rdmvec, rdmvec)
    normrdmvec = arr(rdmvec) / math.sqrt(scale)
    mat.append(normrdmvec)
    return mat


def ortho(dim, matrix, start):
    if dim == 1:
        test_vec = matrix
        scale = dot(test_vec, test_vec)
        normvec = arr(test_vec) / math.sqrt(scale)
        return normvec
    else:
        for i in range(start, dim):
            test_vec = matrix[i]
            scale = dot(test_vec, test_vec)
            normvec = arr(test_vec) / math.sqrt(scale)
            curr_vec = matrix[0]
            scale = dot(curr_vec, curr_vec)
            normcurr_vec = arr(curr_vec) / math.sqrt(scale)
            sum_proj = proj(arr(normcurr_vec), arr(normvec))
            for j in range(1, dim):
                if j == i:
                    continue
                curr_vec = matrix[j]
                scale = dot(curr_vec, curr_vec)
                norm_currvec = arr(curr_vec) / math.sqrt(scale)
                sum_proj += proj(norm_currvec, arr(normvec))
            finalvec = arr(normvec) - arr(sum_proj)
            scale = dot(finalvec, finalvec)
            normfinalvec = arr(finalvec) / math.sqrt(scale)
            matrix[i] = normfinalvec
        return matrix


def checknorm(matrix):
    sumvec = []
    for i in range(0, len(matrix)):
        sum = 0.0
        for j in range(0, len(matrix)):
            norm = dot(matrix[i], matrix[j])
            absnorm = sqr(norm * norm)
            sum += absnorm
        sumvec.append(sum)
    return sumvec


def nma_3Nminus6dof_asfunction(hessfile, basedir):
    extr = imp.load_source(
        "orcaextractors", str(basedir) + "/operations/hessian_from_hes.py"
    )
    extrmass = imp.load_source(
        "orcaextractors", str(basedir) + "/operations/mass_from_hes.py"
    )
    extrxyz = imp.load_source(
        "orcaextractors", str(basedir) + "/operations/geo_from_hes.py"
    )

    hessian = extr.hessian_from_hes(hessfile)
    mass = extrmass.mass_from_hes(hessfile)
    coords = extrxyz.geo_from_hes(hessfile, "B")

    # Begin Center of mass
    com = []
    for i in range(0, 3):
        curr_com = float(0.0)
        sum_mass = float(0.0)
        for j in range(0, len(mass)):
            curr_com += float(mass[j]) * float(coords[j * 3 + i])
            sum_mass += float(mass[j])
        com.append(float(curr_com) / float(sum_mass))
    # End Center of mass
    # Begin shift coords
    shift_coords = coords
    for i in range(0, len(mass)):
        for j in range(0, 3):
            shift_coords[i * 3 + j] = float(shift_coords[i * 3 + j]) - float(com[j])

    t_inert = []
    t_row = []
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele += float(mass[i]) * (
            shift_coords[i * 3 + 1] * shift_coords[i * 3 + 1]
            + shift_coords[i * 3 + 2] * shift_coords[i * 3 + 2]
        )
    t_row.append(ele)
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele -= float(mass[i]) * (shift_coords[i * 3 + 0] * shift_coords[i * 3 + 1])
    t_row.append(ele)
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele -= float(mass[i]) * (shift_coords[i * 3 + 0] * shift_coords[i * 3 + 2])
    t_row.append(ele)
    t_inert.append(t_row)
    t_row = []
    t_row.append(float(t_inert[0][1]))
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele += float(mass[i]) * (
            shift_coords[i * 3 + 0] * shift_coords[i * 3 + 0]
            + shift_coords[i * 3 + 2] * shift_coords[i * 3 + 2]
        )
    t_row.append(ele)
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele -= float(mass[i]) * (shift_coords[i * 3 + 1] * shift_coords[i * 3 + 2])
    t_row.append(ele)
    t_inert.append(t_row)
    t_row = []
    t_row.append(float(t_inert[0][2]))
    t_row.append(float(t_inert[1][2]))
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele += float(mass[i]) * (
            shift_coords[i * 3 + 0] * shift_coords[i * 3 + 0]
            + shift_coords[i * 3 + 1] * shift_coords[i * 3 + 1]
        )
    t_row.append(ele)
    t_inert.append(t_row)
    # generate principal moments and their eigenvectors
    e_princ, princ = LA.eig(arr(t_inert))
    # begin setup transrot
    transrotmat = []
    for i in range(0, 3):
        transvec = []
        for j in range(0, len(mass)):
            for k in range(0, i):
                transvec.append(float(0.0))
            transvec.append(math.sqrt(float(mass[j])))
            for k in range(i + 1, 3):
                transvec.append(float(0.0))
        transrotmat.append(transvec)
    # begin rotvecs
    rotvec = []
    for i in range(0, len(mass)):
        curr_coords = []
        for j in range(0, 3):
            curr_coords.append(shift_coords[3 * i + j])

        py = dot(arr(curr_coords).astype(float), arr(princ[1]).astype(float))
        pz = dot(arr(curr_coords).astype(float), arr(princ[2]).astype(float))
        for j in range(0, 3):
            rotvec.append(
                (py * princ[2][j] - pz * princ[1][j]) * math.sqrt(float(mass[i]))
            )
    transrotmat.append(rotvec)
    rotvec = []
    for i in range(0, len(mass)):
        curr_coords = []
        for j in range(0, 3):
            curr_coords.append(shift_coords[3 * i + j])

        px = dot(arr(curr_coords).astype(float), arr(princ[0]).astype(float))
        pz = dot(arr(curr_coords).astype(float), arr(princ[2]).astype(float))
        for j in range(0, 3):
            rotvec.append(
                (pz * princ[0][j] - px * princ[2][j]) * math.sqrt(float(mass[i]))
            )
    transrotmat.append(rotvec)
    rotvec = []
    for i in range(0, len(mass)):
        curr_coords = []
        for j in range(0, 3):
            curr_coords.append(shift_coords[3 * i + j])

        px = dot(arr(curr_coords).astype(float), arr(princ[0]).astype(float))
        py = dot(arr(curr_coords).astype(float), arr(princ[1]).astype(float))
        for j in range(0, 3):
            rotvec.append(
                (px * princ[1][j] - py * princ[0][j]) * math.sqrt(float(mass[i]))
            )
    transrotmat.append(rotvec)
    # Begin mass weight Hessianxyz
    for i in range(0, len(mass)):
        for j in range(0, 3):
            for k in range(0, len(mass)):
                for l in range(0, 3):
                    hessian[i * 3 + j][k * 3 + l] /= math.sqrt(
                        float(mass[i]) * float(mass[k])
                    )
                    if (
                        hessian[i * 3 + j][k * 3 + l] < 0.00000000000001
                        and hessian[i * 3 + j][k * 3 + l] > 0.0
                    ):
                        hessian[i * 3 + j][k * 3 + l] = float(0.0)
                    if (
                        hessian[i * 3 + j][k * 3 + l] > -0.00000000000001
                        and hessian[i * 3 + j][k * 3 + l] < 0.0
                    ):
                        hessian[i * 3 + j][k * 3 + l] = float(0.0)
    # end mass weight Hessianxyz
    # begin normalize transrot
    norm_trm = []
    for i in range(0, 6):
        scale = dot(transrotmat[i], transrotmat[i])
        norm_vec = arr(transrotmat[i]) / math.sqrt(scale)
        norm_trm.append(arr(norm_vec).astype(float))
    # end normalize transrot

    # Begin Ortho
    # Ortho for rotvecs
    for i in range(3, 6):
        ortho(6, norm_trm, 3)
    for i in range(6, len(mass) * 3):
        norm_trm = newvec(norm_trm, len(mass) * 3, i)
        ortho(i + 1, norm_trm, i)
    enter = 0
    found = False
    for i in range(0, len(mass) * 3):
        currnorm = float(0.0)
        for j in range(0, len(mass) * 3):
            currnorm += dot(norm_trm[i], norm_trm[j])
        if currnorm > 1.01:
            found = True
    while enter <= 500 and found:
        ortho(len(mass) * 3, norm_trm, 6)
        found = False
        for i in range(0, len(mass) * 3):
            currnorm = float(0.0)
            for j in range(0, len(mass) * 3):
                currnorm += dot(norm_trm[i], norm_trm[j])
            if currnorm > 1.01:
                found = True
        enter += 1
    if enter >= 500:
        print "had to orthogonalize more than 500 times! Check your results!"
    # ENd Ortho
    d_mat = []
    for i in range(6, len(mass) * 3):
        d_mat.append(norm_trm[i])
    dh = dot(d_mat, hessian)
    h_int = dot(dh, transpose(d_mat))
    e_val, nm_matrix = LA.eig(arr(h_int))
    mass_sqrtvec_matrix = []
    for i in range(0, len(mass)):
        for p in range(0, 3):
            massline = []
            for j in range(0, len(mass)):
                for k in range(0, 3):
                    if i == j and p == k:
                        massline.extend([float(1.0 / sqr(float(mass[i])))])
                    else:
                        massline.extend([float(0.0)])
            mass_sqrtvec_matrix.append(massline)
    mw_d_mat = dot(arr(mass_sqrtvec_matrix), transpose(d_mat))
    cart_nm = dot(mw_d_mat, nm_matrix)
    # begin sort eigenvalues/vectors
    done = False
    while not done:
        done = True
        for i in range(0, len(e_val) - 1):
            if float(e_val[i + 1]) < float(e_val[i]):
                done = False
                nm_matrix[:, [i, i + 1]] = nm_matrix[:, [i + 1, i]]
                cart_nm[:, [i, i + 1]] = cart_nm[:, [i + 1, i]]
                e_val[i + 1], e_val[i] = e_val[i], e_val[i + 1]
    return e_val, nm_matrix


def nma_3Nminus6dof(sysargs):
    basedir = sysargs[2]
    extr = imp.load_source(
        "orcaextractors", str(basedir) + "/operations/hessian_from_hes.py"
    )
    extrmass = imp.load_source(
        "orcaextractors", str(basedir) + "/operations/mass_from_hes.py"
    )
    extrxyz = imp.load_source(
        "orcaextractors", str(basedir) + "/operations/geo_from_hes.py"
    )

    hessian = extr.hessian_from_hes(sysargs[1])
    mass = extrmass.mass_from_hes(sysargs[1])
    coords = extrxyz.geo_from_hes(sysargs[1], "B")

    # Begin Center of mass
    com = []
    for i in range(0, 3):
        curr_com = float(0.0)
        sum_mass = float(0.0)
        for j in range(0, len(mass)):
            curr_com += float(mass[j]) * float(coords[j * 3 + i])
            sum_mass += float(mass[j])
        com.append(float(curr_com) / float(sum_mass))
    # End Center of mass
    # Begin shift coords
    shift_coords = coords
    for i in range(0, len(mass)):
        for j in range(0, 3):
            shift_coords[i * 3 + j] = float(shift_coords[i * 3 + j]) - float(com[j])
    # coords still match! good!
    # End shift coords
    # begin inertia tensor
    t_inert = []
    t_row = []
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele += float(mass[i]) * (
            shift_coords[i * 3 + 1] * shift_coords[i * 3 + 1]
            + shift_coords[i * 3 + 2] * shift_coords[i * 3 + 2]
        )
    t_row.append(ele)
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele -= float(mass[i]) * (shift_coords[i * 3 + 0] * shift_coords[i * 3 + 1])
    t_row.append(ele)
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele -= float(mass[i]) * (shift_coords[i * 3 + 0] * shift_coords[i * 3 + 2])
    t_row.append(ele)
    t_inert.append(t_row)
    t_row = []
    t_row.append(float(t_inert[0][1]))
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele += float(mass[i]) * (
            shift_coords[i * 3 + 0] * shift_coords[i * 3 + 0]
            + shift_coords[i * 3 + 2] * shift_coords[i * 3 + 2]
        )
    t_row.append(ele)
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele -= float(mass[i]) * (shift_coords[i * 3 + 1] * shift_coords[i * 3 + 2])
    t_row.append(ele)
    t_inert.append(t_row)
    t_row = []
    t_row.append(float(t_inert[0][2]))
    t_row.append(float(t_inert[1][2]))
    ele = float(0.0)
    for i in range(0, len(mass)):
        ele += float(mass[i]) * (
            shift_coords[i * 3 + 0] * shift_coords[i * 3 + 0]
            + shift_coords[i * 3 + 1] * shift_coords[i * 3 + 1]
        )
    t_row.append(ele)
    t_inert.append(t_row)
    # end inertia tensor
    # generate principal moments and their eigenvectors
    e_princ, princ = LA.eig(arr(t_inert))
    # e_princ ok, princ ok!
    # begin setup transrot
    transrotmat = []
    for i in range(0, 3):
        transvec = []
        for j in range(0, len(mass)):
            for k in range(0, i):
                transvec.append(float(0.0))
            transvec.append(math.sqrt(float(mass[j])))
            for k in range(i + 1, 3):
                transvec.append(float(0.0))
        transrotmat.append(transvec)
    # begin rotvecs
    rotvec = []
    for i in range(0, len(mass)):
        curr_coords = []
        for j in range(0, 3):
            curr_coords.append(shift_coords[3 * i + j])
        py = dot(arr(curr_coords).astype(float), arr(princ[1]).astype(float))
        pz = dot(arr(curr_coords).astype(float), arr(princ[2]).astype(float))
        for j in range(0, 3):
            rotvec.append(
                (py * princ[2][j] - pz * princ[1][j]) * math.sqrt(float(mass[i]))
            )
    transrotmat.append(rotvec)
    rotvec = []
    for i in range(0, len(mass)):
        curr_coords = []
        for j in range(0, 3):
            curr_coords.append(shift_coords[3 * i + j])
        # 	print curr_coords
        px = dot(arr(curr_coords).astype(float), arr(princ[0]).astype(float))
        pz = dot(arr(curr_coords).astype(float), arr(princ[2]).astype(float))
        for j in range(0, 3):
            rotvec.append(
                (pz * princ[0][j] - px * princ[2][j]) * math.sqrt(float(mass[i]))
            )
    transrotmat.append(rotvec)
    rotvec = []
    for i in range(0, len(mass)):
        curr_coords = []
        for j in range(0, 3):
            curr_coords.append(shift_coords[3 * i + j])
        # 	print curr_coords
        px = dot(arr(curr_coords).astype(float), arr(princ[0]).astype(float))
        py = dot(arr(curr_coords).astype(float), arr(princ[1]).astype(float))
        for j in range(0, 3):
            rotvec.append(
                (px * princ[1][j] - py * princ[0][j]) * math.sqrt(float(mass[i]))
            )
    transrotmat.append(rotvec)
    # Begin mass weight Hessianxyz
    for i in range(0, len(mass)):
        for j in range(0, 3):
            for k in range(0, len(mass)):
                for l in range(0, 3):
                    hessian[i * 3 + j][k * 3 + l] /= math.sqrt(
                        float(mass[i]) * float(mass[k])
                    )
                    if (
                        hessian[i * 3 + j][k * 3 + l] < 0.00000000000001
                        and hessian[i * 3 + j][k * 3 + l] > 0.0
                    ):
                        hessian[i * 3 + j][k * 3 + l] = float(0.0)
                    if (
                        hessian[i * 3 + j][k * 3 + l] > -0.00000000000001
                        and hessian[i * 3 + j][k * 3 + l] < 0.0
                    ):
                        hessian[i * 3 + j][k * 3 + l] = float(0.0)
    # end mass weight Hessianxyz
    # begin normalize transrot
    norm_trm = []
    for i in range(0, 6):
        scale = dot(transrotmat[i], transrotmat[i])
        norm_vec = arr(transrotmat[i]) / math.sqrt(scale)
        norm_trm.append(arr(norm_vec).astype(float))
    # end normalize transrot
    # Begin Ortho
    # Ortho for rotvecs
    for i in range(3, 6):
        ortho(6, norm_trm, 3)
    for i in range(6, len(mass) * 3):
        norm_trm = newvec(norm_trm, len(mass) * 3, i)
        ortho(i + 1, norm_trm, i)
    enter = 0
    found = False
    for i in range(0, len(mass) * 3):
        currnorm = float(0.0)
        for j in range(0, len(mass) * 3):
            currnorm += dot(norm_trm[i], norm_trm[j])
        if currnorm > 1.01:
            found = True
    while (enter <= 500 and found) or enter < 10:
        ortho(len(mass) * 3, norm_trm, 6)
        found = False
        for i in range(0, len(mass) * 3):
            currnorm = float(0.0)
            for j in range(0, len(mass) * 3):
                currnorm += dot(norm_trm[i], norm_trm[j])
            if 1.0 - (sqr(currnorm * currnorm)) > 0.01:
                found = True
        enter += 1
    if enter >= 500:
        print "had to orthogonalize more than 500 times! Check your results!"
    # ENd Ortho
    d_mat = []
    for i in range(6, len(mass) * 3):
        d_mat.append(norm_trm[i])
    dh = dot(d_mat, hessian)
    h_int = dot(dh, transpose(d_mat))
    e_val, nm_matrix = LA.eigh(arr(hessian))
    b_mat = construct_bmat(shift_coords, basedir)
    # begin sort eigenvalues/vectors
    done = False
    while not done:
        done = True
        for i in range(0, len(e_val) - 1):
            if float(e_val[i + 1]) < float(e_val[i]):
                done = False
                nm_matrix[:, [i, i + 1]] = nm_matrix[:, [i + 1, i]]
                e_val[i + 1], e_val[i] = e_val[i], e_val[i + 1]
    for i in range(0, len(e_val)):
        if e_val[i] < 0.0:
            print str(-1.0 * sqr(abs(float(e_val[i])) * float(26424608)))
        else:
            print str(sqr(float(e_val[i]) * float(26424608)))
    exit(1)
    print h_int
    print nm_matrix
    print arr(nm_matrix)


if __name__ == "__main__":
    import sys

    nma_3Nminus6dof(sys.argv)
