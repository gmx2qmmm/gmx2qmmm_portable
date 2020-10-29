#!/usr/bin/env python2
# encoding: ISO-8859-15

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# This will take a Turbomole point_charges file and generate sum of charge vectors at a given point in space

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$"  # during a rain storm

import re
import sys

import numpy as np


def sum_pcf_tm(inp, x, y, z):
    base = [float(x), float(y), float(z)]
    sumvec = np.array([0.0, 0.0, 0.0])
    with open(inp) as ifile:
        for line in ifile:
            match = re.search(
                r"^\s*([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                new_vec = [
                    float(match.group(1)),
                    float(match.group(2)),
                    float(match.group(3)),
                ]
                new_charge = float(match.group(4))
                distvec = np.array(new_vec) - np.array(base)
                qdistvec = np.array(distvec) * new_charge
                sumvec += qdistvec
    return sumvec


def sum_pcf_tm_nofile(inp, x, y, z):
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


def sum_pcf_tm_red(inp, x, y, z, atomlist):
    base = [float(x), float(y), float(z)]
    sumvec = np.array([0.0, 0.0, 0.0])
    count = 0
    with open(inp) as ifile:
        for line in ifile:
            count += 1
            if count not in atomlist:
                continue
            match = re.search(
                r"^\s*([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                new_vec = [
                    float(match.group(1)),
                    float(match.group(2)),
                    float(match.group(3)),
                ]
                new_charge = float(match.group(4))
                if count != 5926:
                    new_charge += 0.233333333
                distvec = np.array(new_vec) - np.array(base)
                qdistvec = np.array(distvec) * new_charge
                sumvec += qdistvec
    return sumvec


if __name__ == "__main__":
    sumvec = sum_pcf_tm(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print(sumvec)
