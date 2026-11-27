import numpy as np
from pytest import approx
import pathlib 

from gmx2qmmm.generators.geometry import * 
DATA = pathlib.Path(__file__).parent / "data"

# Testing data. Coordinates in geom.g96/geom.gro are in nm.
# geom.g96 contains 3 not connected atoms (Gly N, Ser CA, Ser CB).
expected_geometry = [47.83, 85.62, 31.12, 49.14, 88.15, 34.21, 50.38, 88.90, 33.76]
input_coords =  [0.123, 0.456, 0.789, 0.321, 0.654, 0.987]
input_charges = [-0.57, +1.57]

# utility function
def compare_lists_of_floats(list_a, list_b):
    return all([approx(a) == b for a, b in zip(list_a, list_b)])


def test_readg96():
    """Tests reading .g96 files and converting units to Angstrom."""
    coords = readg96(DATA / "geom.g96")
    assert len(coords) == len(expected_geometry)
    assert compare_lists_of_floats(coords, expected_geometry)


def test_readgeo():
    """Tests reading .gro files and converting to Angstrom."""
    coords = readgeo(DATA / "geom.gro")
    assert len(coords) == len(expected_geometry)
    assert compare_lists_of_floats(coords, expected_geometry)


def test_make_xyzq():
    """Tests building an Nx4 array [x, y, z, q] from a flat coordinate list and charges."""
    res_xyzq = make_xyzq(input_coords, input_charges)
    ref_xyzq = np.array([
        [0.123, 0.456, 0.789, -0.57],
        [0.321, 0.654, 0.987, +1.57],
    ])
    assert np.allclose(res_xyzq, ref_xyzq)


def test_make_xyzq_io():
    """Tests that charges of outer-list atoms (1-based) are set to 0."""
    res_xyzq = make_xyzq_io(input_coords, input_charges, [1])
    ref_xyzq = np.array([
        [0.123, 0.456, 0.789, 0.0],
        [0.321, 0.654, 0.987, +1.57],
    ])
    assert np.allclose(res_xyzq, ref_xyzq)


def test_read_gmx_structure_atoms():
    """Tests reading atom labels from .g96 files."""
    atom_info = read_gmx_structure_atoms(DATA / "geom.g96")
    atom_info_ref = ["   75 GLY       N      1", "   76 SER      CA     12", "   76 SER      CB     14"]
    assert len(atom_info) == len(atom_info_ref)
    assert atom_info == atom_info_ref 


def test_read_gmx_box_vectors():
    """Tests reading box vectors in nm."""
    box_info = read_gmx_box_vectors(DATA / "geom.g96")
    box_info_ref = [21.14, 12.8, 6.069]
    assert len(box_info) == len(box_info_ref)
    assert compare_lists_of_floats(box_info, box_info_ref)


def test_read_gmx_structure_header():
    """Tests reading the first four header lines verbatim."""
    header = read_gmx_structure_header(DATA / "geom.g96") 
    header_ref = 'TITLE\nProtein\nEND\nPOSITION\n'
    assert header == header_ref


def test_write_g96(tmp_path):
    """Tests writing .g96 files (Angstrom to nm) against the reference file."""
    ref_file = DATA / "geom.g96"
    new_file = tmp_path / "test.g96"
    write_g96( new_file 
             , read_gmx_structure_header(ref_file)
             , read_gmx_structure_atoms(ref_file)
             , np.array(expected_geometry).reshape(-1, 3) 
             , read_gmx_box_vectors(ref_file)
             )

    # compare line by line omitting "\n"
    with open(new_file, 'r') as a, open(ref_file, 'r') as b:
        for line_a, line_b in zip(a,b):
            assert line_a.rstrip() == line_b.rstrip()


def test_write_and_read_agree(tmp_path):
    """Tests that coordinates and box survive a write_g96/readg96 round trip."""
    tmp_file = tmp_path / "tmp.g96"
    data_file = DATA / "geom.g96"
    write_g96(tmp_file,
              read_gmx_structure_header(data_file),
              read_gmx_structure_atoms(data_file),
              np.array(expected_geometry).reshape(-1, 3),
              read_gmx_box_vectors(data_file))
    assert compare_lists_of_floats(readg96(tmp_file), expected_geometry)
    assert compare_lists_of_floats(read_gmx_box_vectors(tmp_file), [21.14, 12.8, 6.069])


def test_readg96_and_readgeo_agree():
    """Tests that .g96 and .gro readers give the same coordinates."""
    assert compare_lists_of_floats(readg96(DATA / "geom.g96"), readgeo(DATA / "geom.gro"))


def test_propagate_dispvec_steep_being_first():
    """Tests that any propagator falls back to steepest descent at curr_step == 1."""
    forces = np.array([[0.0, 0.0, 2.0]])
    steep = propagate_dispvec("STEEP", None, [forces], 2.0, 1.0, curr_step=1)
    cg    = propagate_dispvec("CONJGRAD", None, [forces], 2.0, 1.0, curr_step=1)
    assert np.allclose(steep, cg)


def test_propagate_dispvec_steep():
    """Tests steepest descent: largest displacement equals stepsize (Bohr to Angstrom)."""
    forces = np.array([[0.0, 0.0, 2.0], [0.0, -1.0, 0.0]])
    disp = propagate_dispvec("STEEP", None, [forces], -2.0, 1.0, curr_step=5)
    bohr = 0.52917721
    assert np.allclose(disp, [[0, 0, bohr], [0, -bohr / 2, 0]])



