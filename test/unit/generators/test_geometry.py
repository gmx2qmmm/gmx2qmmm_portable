import numpy as np
from pytest import approx
import pathlib 

from gmx2qmmm.generators.geometry import * 

# Testing data
expected_geometry = [47.83, 85.62, 31.12, 49.14, 88.15, 34.21, 50.38, 88.90, 33.76]
input_coords =  [0.123, 0.456, 0.789, 0.321, 0.654, 0.987]
input_charges = [-0.57, +1.57]

# utility function
def compare_lists_of_floats(list_a, list_b):
    return all([approx(a) == b for a, b in zip(list_a, list_b)])


def test_readg96():
    coords = readg96("test/unit/data/geom.g96")
    assert len(coords) == len(expected_geometry)
    assert compare_lists_of_floats(coords, expected_geometry)


def test_readgeo():
    coords = readgeo("test/unit/data/geom.gro")
    assert len(coords) == len(expected_geometry)
    assert compare_lists_of_floats(coords, expected_geometry)


def test_make_xyzq():
    res_xyzq = make_xyzq(input_coords, input_charges)
    ref_xyzq = np.array([
        [0.123, 0.456, 0.789, -0.57],
        [0.321, 0.654, 0.987, +1.57],
    ])
    assert np.allclose(res_xyzq, ref_xyzq)


def test_make_xyzq_io():
    res_xyzq = make_xyzq_io(input_coords, input_charges, [1])
    ref_xyzq = np.array([
        [0.123, 0.456, 0.789, 0.0],
        [0.321, 0.654, 0.987, +1.57],
    ])
    assert np.allclose(res_xyzq, ref_xyzq)


def test_read_gmx_structure_atoms():
    atom_info = read_gmx_structure_atoms("test/unit/data/geom.g96")
    atom_info_ref = ["   75 GLY       N      1", "   76 SER      CA     12", "   76 SER      CB     14"]
    assert len(atom_info) == len(atom_info_ref)
    assert atom_info == atom_info_ref 


def test_read_gmx_box_vectors():
    box_info = read_gmx_box_vectors("test/unit/data/geom.g96")
    box_info_ref = [21.14, 12.8, 6.069]
    assert len(box_info) == len(box_info_ref)
    assert compare_lists_of_floats(box_info, box_info_ref)


def test_read_gmx_structure_header():
    header = read_gmx_structure_header("test/unit/data/geom.g96") 
    header_ref = 'TITLE\nProtein\nEND\nPOSITION\n'
    assert header == header_ref


def test_write_g96():
    ref_file = "test/unit/data/geom.g96"
    new_file = "test/unit/data/test.g96"
    write_g96( new_file 
             , read_gmx_structure_header(ref_file)
             , read_gmx_structure_atoms(ref_file)
             , [[47.83, 85.62, 31.12], [49.14, 88.15, 34.21], [50.38, 88.90, 33.76]] 
             , read_gmx_box_vectors(ref_file)
             )

    # compare line by line omitting "\n"
    with open(new_file, 'r') as a, open(ref_file, 'r') as b:
        for line_a, line_b in zip(a,b):
            assert line_a.rstrip() == line_b.rstrip()
    
    # remove temporary file
    file_to_rem = pathlib.Path(new_file)
    file_to_rem.unlink()

