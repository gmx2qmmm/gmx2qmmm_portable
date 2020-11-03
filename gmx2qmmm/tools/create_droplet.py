#!/usr/bin/env python2
# encoding: ISO-8859-15

# This will take a cubic simulation box and create a water droplet with a given radius for QM/MM simulation using gmx2qmmm.
# NEEDS a cubic box. Will maintain charge of original system by keeping all elements of the system except for solvent.
# Attention: Ions larger than one atom will be considered a non-ion residue. This might lead to error messages. but the program is designed to print the coordinates anyway.
# Needs input .gro file, input .top file, radius (angstrom), number of atom to center the droplet on, name of the solvent residue (i.e., SOL or H2O), output .gro file, output .top file.


def find_center_coords(gro, n_center):
    import re

    coords = []
    with open(gro) as ifile:
        counter = -2
        for line in ifile:
            counter += 1
            if counter == int(n_center):
                match = re.search(
                    r"^\s*\d+\S+\s+\S+\s*\d+\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    print "Something is wrong with your .gro file. Maybe the format is unexpected? Exiting. Last line:"
                    print line
                    exit(1)
                coords = [
                    float(match.group(1)) * 10.0,
                    float(match.group(2)) * 10.0,
                    float(match.group(3)) * 10.0,
                ]
                break
    return coords


def analyze_gro(gro):
    import re

    coords = []
    residues = []
    box_vectors = [0.0, 0.0, 0.0]
    # names=[]
    n_a = 0
    with open(gro) as ifile:
        counter = -2
        for line in ifile:
            counter += 1
            if counter == 0:
                match = re.search(r"^\s*(\d+)", line, flags=re.MULTILINE)
                if not match:
                    print "Something is wrong with your .gro file. Maybe the format is unexpected? Exiting. Last line:"
                    print line
                    exit(1)
                n_a = int(match.group(1))
            if counter > 0 and counter <= n_a:
                match = re.search(
                    r"^(.{5})(.{5})(.{5})\s*\d+\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    print "Something is wrong with your .gro file. Maybe the format is unexpected? Exiting. Last line:"
                    print line
                    exit(1)
                resline = [
                    int(match.group(1)),
                    (match.group(2)).strip(),
                    (match.group(3)).strip(),
                ]
                residues.append(resline)
                coordline = [
                    float(match.group(4)) * 10.0,
                    float(match.group(5)) * 10.0,
                    float(match.group(6)) * 10.0,
                ]
                coords.append(coordline)
            if counter > n_a:
                match = re.search(
                    r"^\s*(\d+\.\d*)\s*(\d+\.\d*)\s*(\d+\.\d*)",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    print "Something is wrong with your .gro file. Maybe the format is unexpected? Exiting. Last line:"
                    print line
                    exit(1)
                box_vectors = [
                    float(match.group(1)) * 10.0,
                    float(match.group(2)) * 10.0,
                    float(match.group(3)) * 10.0,
                ]

    if len(coords) != n_a:
        print str(
            "Expected "
            + str(n_a)
            + " atoms. Found "
            + str(len(n_a))
            + " atoms. Exiting."
        )
        exit(1)
    if (box_vectors[0] + box_vectors[1] + box_vectors[2]) < 0.00001:
        print str("Strange box encountered. Exiting. Box vectors were:")
        print box_vectors
        exit(1)
    return coords, residues, box_vectors


def recenter_coords(coords, center):
    new_coords = []
    for element in coords:
        new_coordline = []
        for i in range(0, 3):
            new_coordline.append(float(float(element[i]) - float(center[i])))
        new_coords.append(new_coordline)
    return new_coords


def check_max_r_for_all_mirror_images(coordline, box_vectors, max_r):
    from numpy import array as arr
    from numpy import linalg as LA

    curr_max = float(max_r)
    # print coordline
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_max > LA.norm(curr_coords):
        curr_max = LA.norm(curr_coords)
    return float(curr_max)


def get_smallest_distance_to_mirror_images(coords, residues, solname, box_vectors):
    max_r = 0.0
    for element in box_vectors:
        if float(element) > max_r:
            max_r = float(element)
    for i in range(0, len(coords) - 1):
        if residues[i][1] == solname:
            continue
        if residues[i][0] == residues[i + 1][0]:

            max_r = check_max_r_for_all_mirror_images(coords[i], box_vectors, max_r)
        if i > 0 and residues[i][0] != residues[i + 1][0]:
            if residues[i][0] == residues[i - 1][0]:
                max_r = check_max_r_for_all_mirror_images(coords[i], box_vectors, max_r)
    if residues[len(coords) - 1][0] == residues[len(coords) - 2][0]:
        max_r = check_max_r_for_all_mirror_images(
            coords[len(coords) - 1], box_vectors, max_r
        )
    return max_r


def find_max_dropsize(gro, r, center, solname):

    coords, residues, box_vectors = analyze_gro(gro)

    coords = recenter_coords(coords, center)
    max_r = get_smallest_distance_to_mirror_images(
        coords, residues, solname, box_vectors
    )
    if max_r < float(r):
        print "Your requested drop radius intersects with mirror images of a large non-solvent compound."
        print "Consequently, your droplet is not a droplet. Will print the coordinates anyway, but consider changing your drop size or recreating your MD with a larger box size."

    return coords, residues, box_vectors


def verify_solvent_coords(coordline, ra, box_vectors):
    from numpy import array as arr
    from numpy import linalg as LA

    r = float(ra)
    disp = []
    curr_coords = arr(coordline)
    if r >= LA.norm(curr_coords):
        disp.append([0, 0, 0])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, 0, 0])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, 0, 0])
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([0, 1, 0])
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([0, -1, 0])
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([0, 0, 1])
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([0, 0, -1])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, 1, 0])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, 1, 0])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, -1, 0])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, -1, 0])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, 0, 1])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, 0, 1])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, 0, -1])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, 0, -1])
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([0, 1, 1])
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([0, -1, 1])
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([0, 1, -1])
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([0, -1, -1])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, 1, 1])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, 1, 1])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, -1, 1])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, 1, -1])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, -1, 1])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, 1, -1])
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([1, -1, -1])
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if r >= LA.norm(curr_coords):
        disp.append([-1, -1, -1])
    return disp


def make_solvent_coords(solcoords, disp_vec, box_vectors):
    from numpy import array as arr
    from numpy import linalg as LA

    new_coords = []
    for element in disp_vec:
        for stuff in solcoords:
            new_coords.append(
                arr(element).astype(float) * arr(box_vectors).astype(float)
                + arr(stuff).astype(float)
            )
    return new_coords


def check_min_r_for_all_mirror_images(coordline, box_vectors):
    from numpy import array as arr
    from numpy import linalg as LA

    # print coordline
    curr_coords = arr(coordline)
    curr_min = LA.norm(curr_coords)
    min_coords = coordline
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            0.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            0.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            0.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    curr_coords = arr(coordline) + arr(
        [
            -1.0 * float(box_vectors[0]),
            -1.0 * float(box_vectors[1]),
            -1.0 * float(box_vectors[2]),
        ]
    )
    if curr_min > LA.norm(curr_coords):
        curr_min = LA.norm(curr_coords)
        min_coords = curr_coords
    return min_coords


def build_drop_gro(r, centered_coords, residues, box_vectors, solname):
    from compiler.ast import flatten

    drop_res = []
    drop_coords = []
    new_boxvectors = [0.0, 0.0, 0.0]
    curr_res = int(residues[0][0])
    i = 0
    while i < len(centered_coords) - 1:
        new_coords = []
        new_res = []
        if residues[i][1] == solname:
            solcoords = []
            solcoords.append(centered_coords[i])
            disp_vec = verify_solvent_coords(centered_coords[i], r, box_vectors)
            while residues[i][0] == residues[i + 1][0] and i < len(centered_coords) - 1:
                i += 1
                solcoords.append(centered_coords[i])
                disp = verify_solvent_coords(centered_coords[i], r, box_vectors)
                for element in disp:
                    if element not in disp_vec:
                        disp_vec.append(element)
            new_coords = make_solvent_coords(solcoords, disp_vec, box_vectors)
            for j in range(0, len(disp_vec)):
                for k in range(0, len(solcoords)):
                    new_res.append(
                        [
                            int(curr_res) + j,
                            solname,
                            residues[i + k - len(new_coords) + 1][2],
                        ]
                    )
            for j in range(0, len(new_coords)):
                drop_coords.append(new_coords[j])
                drop_res.append(new_res[j])
            curr_res = (drop_res[len(drop_res) - 1][0]) + 1
            i += 1
            continue
        if residues[i][0] == residues[i + 1][0]:
            new_coords = centered_coords[i]
            new_res = [curr_res, residues[i][1], residues[i][2]]
        elif i > 0 and residues[i][0] != residues[i + 1][0]:
            if residues[i][0] == residues[i - 1][0]:
                new_coords = centered_coords[i]
                new_res = [curr_res, residues[i][1], residues[i][2]]
                curr_res += 1
            else:
                new_coords = centered_coords[i]
                new_res = [curr_res, residues[i][1], residues[i][2]]
                curr_res += 1
        else:  # is NOT a large residue
            new_coords = check_min_r_for_all_mirror_images(
                centered_coords[i], box_vectors
            )
            new_res = [curr_res, residues[i][1], residues[i][2]]
            curr_res += 1
        if len(flatten(new_coords)) == 3:
            drop_coords.append(new_coords)
            drop_res.append(new_res)
        else:
            drop_coords.extend(new_coords)
            drop_res.extend(new_res)
        i += 1
    if residues[len(centered_coords) - 1][1] != solname:
        if (
            residues[len(centered_coords) - 1][0]
            == residues[len(centered_coords) - 2][0]
        ):
            drop_coords.append(centered_coords[len(centered_coords) - 1])
            drop_res.append(
                [
                    curr_res,
                    residues[len(residues) - 1][1],
                    residues[len(residues) - 1][2],
                ]
            )
            curr_res += 1
        else:  # is NOT a large residue
            new_coords = check_min_r_for_all_mirror_images(
                centered_coords[len(centered_coords) - 1], box_vectors
            )
            drop_coords.append(new_coords)
            drop_res.append(
                [
                    curr_res,
                    residues[len(residues) - 1][1],
                    residues[len(residues) - 1][2],
                ]
            )
            curr_res += 1
    for i in range(0, 3):
        min_dist = 0.0
        max_dist = 0.0
        for j in range(0, len(drop_coords)):
            if drop_coords[j][i] < min_dist:
                min_dist = drop_coords[j][i]
            if drop_coords[j][i] > max_dist:
                max_dist = drop_coords[j][i]
        new_boxvectors[i] = float(max_dist) - float(min_dist)
    return drop_coords, drop_res, new_boxvectors


def write_drop_gro(gro, new_coords, new_res, new_vectors, solname, groout):
    import re

    with open(groout, "w") as ofile:
        with open(gro) as ifile:
            counter = -1
            for line in ifile:
                ofile.write(line)
                break
            for line in ifile:
                ofile.write(str(len(new_coords)) + "\n")
                break
            for line in ifile:
                match = re.search(
                    r"^\s*(\d+\.\d*)\s*(\d+\.\d*)\s*(\d+\.\d*)",
                    line,
                    flags=re.MULTILINE,
                )
                if match:
                    ofile.write(
                        "{:>10.5f}{:>10.5f}{:>10.5f}".format(
                            float(new_vectors[0]) / 10.0,
                            float(new_vectors[1]) / 10.0,
                            float(new_vectors[2]) / 10.0,
                        )
                    )
                    continue
                counter += 1
                if counter == len(new_coords):
                    break
                match = re.search(
                    r"^(.{5})(.{5})(.{5})(.{5})\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})",
                    line,
                    flags=re.MULTILINE,
                )
                if match:
                    if (match.group(2)).strip() != solname:
                        index2 = counter + 1 % 100000
                        ofile.write(
                            "{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(
                                int(new_res[counter][0]) % 100000,
                                match.group(2),
                                match.group(3),
                                index2,
                                new_coords[counter][0] / 10.0,
                                new_coords[counter][1] / 10.0,
                                new_coords[counter][2] / 10.0,
                            )
                        )
                        continue
                    else:
                        for line in ifile:
                            match = re.search(
                                r"^\s*(\d+\.\d*)\s*(\d+\.\d*)\s*(\d+\.\d*)",
                                line,
                                flags=re.MULTILINE,
                            )
                            if match:
                                ofile.write(
                                    "{:>10.5f}{:>10.5f}{:>10.5f}".format(
                                        float(new_vectors[0]) / 10.0,
                                        float(new_vectors[1]) / 10.0,
                                        float(new_vectors[2]) / 10.0,
                                    )
                                )
                                break
                            match = re.search(
                                r"^(.{5})(.{5})(.{5})(.{5})\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})",
                                line,
                                flags=re.MULTILINE,
                            )
                            if not match:
                                print "Error reading gro file. Last line:"
                                print line
                                exit(1)
                            if (match.group(2)).strip() != solname:
                                break
                        while new_res[counter][1] == solname:
                            # index=int(match.group(1))%100000
                            index2 = counter + 1 % 100000
                            # print str(counter+1)
                            # print new_coords[counter]
                            # print new_res[counter]
                            ofile.write(
                                "{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(
                                    int(new_res[counter][0]) % 100000,
                                    new_res[counter][1],
                                    new_res[counter][2],
                                    index2,
                                    new_coords[counter][0] / 10.0,
                                    new_coords[counter][1] / 10.0,
                                    new_coords[counter][2] / 10.0,
                                )
                            )
                            counter += 1
                            if counter == len(new_res):
                                counter -= 1
                                break
                        match = re.search(
                            r"^\s*(\d+\.\d*)\s*(\d+\.\d*)\s*(\d+\.\d*)",
                            line,
                            flags=re.MULTILINE,
                        )
                        if match:
                            ofile.write(
                                "{:>10.5f}{:>10.5f}{:>10.5f}".format(
                                    float(new_vectors[0]) / 10.0,
                                    float(new_vectors[1]) / 10.0,
                                    float(new_vectors[2]) / 10.0,
                                )
                            )
                            continue
                        match = re.search(
                            r"^(.{5})(.{5})(.{5})(.{5})\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})\s*([-]*\d+\.\d{3})",
                            line,
                            flags=re.MULTILINE,
                        )
                        if match:
                            if (match.group(2)).strip() != solname:
                                # counter+=1
                                index2 = counter + 1 % 100000
                                ofile.write(
                                    "{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(
                                        int(new_res[counter][0]) % 100000,
                                        match.group(2),
                                        match.group(3),
                                        index2,
                                        new_coords[counter][0] / 10.0,
                                        new_coords[counter][1] / 10.0,
                                        new_coords[counter][2] / 10.0,
                                    )
                                )
                        continue
                print "Unexpected line in input gro file. Exiting. Last line:"
                print line
                exit(1)
        ofile.write("\n")


def count_sol(new_res, solname):
    n_sol = 0
    old_res = 0
    for element in new_res:
        curr_res = element[0]
        if curr_res == old_res:
            continue
        if element[1].strip() == solname:
            n_sol += 1
            old_res = curr_res
    return n_sol


def update_mollist(mollist, solname, n_sol):
    new_mollist = []
    found = False
    for element in mollist:
        if element[0].strip() == solname and found:
            continue
        new_mollist.append(element)
        if element[0].strip() == solname and not found:
            new_mollist[len(new_mollist) - 1][1] = n_sol
            found = True
    return new_mollist


def write_drop_top(top, new_mollist, topout):
    import re

    with open(topout, "w") as ofile:
        with open(top) as ifile:
            found1 = False
            for line in ifile:
                match = re.search(r"^\s*\[\s*molecules\s*\]", line, flags=re.MULTILINE)
                if not match:
                    ofile.write(line)
                else:
                    if found1:
                        print "Found two molecules entries in top file. Please clean up first. Exiting."
                        exit(1)
                    found1 = True
                    ofile.write(line)
                    found = False
                    for line in ifile:
                        match = re.search(r"^\s*;", line, flags=re.MULTILINE)
                        if match:
                            ofile.write(line)
                            continue
                        match = re.search(r"^\s*\n", line, flags=re.MULTILINE)
                        if match:
                            break
                        if not match and not found:
                            found = True
                            for element in new_mollist:
                                ofile.write(
                                    "{:<20s}{:>10d}\n".format(
                                        element[0], int(element[1])
                                    )
                                )
                    ofile.write("\n")


def create_droplet(gro, top, r, n_center, solname, groout, topout):
    import imp
    import os.path

    basedir = os.path.dirname(os.path.abspath(__file__))
    make_pcf = imp.load_source(
        "operations", str(basedir + "/../pointcharges/generate_pcf_from_top.py")
    )
    from numpy import array as arr

    center = find_center_coords(gro, n_center)
    centered_coords, residues, box_vectors = find_max_dropsize(gro, r, center, solname)
    new_coords, new_res, new_vectors = build_drop_gro(
        r, centered_coords, residues, box_vectors, solname
    )
    write_drop_gro(gro, new_coords, new_res, new_vectors, solname, groout)
    mollist = make_pcf.readmols(top)
    n_sol = count_sol(new_res, solname)
    new_mollist = update_mollist(mollist, solname, n_sol)
    write_drop_top(top, new_mollist, topout)


if __name__ == "__main__":
    import sys

    create_droplet(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        sys.argv[4],
        sys.argv[5],
        sys.argv[6],
        sys.argv[7],
    )
