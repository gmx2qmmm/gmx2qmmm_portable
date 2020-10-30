import collections
import datetime

import numpy as np


def _flatten(x):
    """Replace deprecated ``compiler.ast.flatten``"""
    for e in x:
        if not isinstance(e, collections.abc.Iterable) or isinstance(e, str):
            yield e
        else:
            yield from _flatten(e)


def logger(log, logstring):
    with open(log, "a") as ofile:
        ofile.write(str(datetime.datetime.now()) + " " + logstring)


def stepper(filename, step):
    """Move to more appropriate module"""
    if step == 0:
        return filename
    else:
        if len(filename) > 7:
            if filename[-7:] == ".fort.7":
                new_filename = str(filename[:-7] + "." + str(step) + ".fort.7")
                return new_filename
        for i in range(0, len(filename)):
            buffer = 0
            if filename[i] == ".":
                buffer = i
        new_filename = str(filename[:buffer] + "." + str(step) + filename[buffer:])
        return new_filename


def make_xyzq(geo, chargevec):
    """Move to more appropriate module"""
    xyzq = []
    count = 0
    for element in chargevec:
        xyzq.append(
            [
                float(geo[count * 3 + 0]),
                float(geo[count * 3 + 1]),
                float(geo[count * 3 + 2]),
                float(element),
            ]
        )
        count += 1
    return xyzq


def get_linkatoms_ang(xyzq, qmatomlist, m1list, connlist, prev_scalefacs):
    """Move to more appropriate module"""
    linkatoms = []
    pairlist = []
    for entry in m1list:
        pair = []
        for i in range(0, len(connlist)):
            if int(entry) in np.array(connlist[i]).astype(int):
                if int(entry) == int(connlist[i][0]):
                    for j in range(0, len(connlist[i])):
                        if int(connlist[i][j]) in np.array(qmatomlist).astype(int):
                            pair = [int(connlist[i][j]), int(entry)]
                            pairlist.append(pair)
                            break
                    break
                else:
                    if int(connlist[i][0]) in np.array(qmatomlist).astype(int):
                        pair = [int(connlist[i][0]), int(entry)]
                        pairlist.append(pair)
                        break
    count = 0
    for entry in pairlist:
        linkcoords = []
        q = [
            float(xyzq[entry[0] - 1][0]),
            float(xyzq[entry[0] - 1][1]),
            float(xyzq[entry[0] - 1][2]),
        ]
        m = [
            float(xyzq[entry[1] - 1][0]),
            float(xyzq[entry[1] - 1][1]),
            float(xyzq[entry[1] - 1][2]),
        ]
        qmvec = np.array(m) - np.array(q)
        scalefac = 0.71290813568205  # ratio between B3LYP/6-31G* optimum of C-C in butane vs C-Link relaxed butane

        if len(prev_scalefacs) != 0:
            scalefac = prev_scalefacs[count][3]
        qmvec_norm_scale = np.array(qmvec) * scalefac
        linkcoords = [
            float(qmvec_norm_scale[0]) + float(q[0]),
            float(qmvec_norm_scale[1]) + float(q[1]),
            float(qmvec_norm_scale[2]) + float(q[2]),
            float(scalefac),
        ]
        linkatoms.append(linkcoords)
        count += 1
    return linkatoms