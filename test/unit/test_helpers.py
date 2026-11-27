import numpy as np
import pytest

from gmx2qmmm.generators._helper import (
    _flatten, filter_xyzq, mask_atoms, normalized_vector,
    get_coordinates_linkatoms_angstrom,
)
from gmx2qmmm.generators.system import SystemInfo

CHAIN_CONN = [[1, 2], [2, 3], [3, 4], [4, 5]]
CHAIN_XYZQ = np.array([[0.0, 0, 0, -0.1],
                       [1.5, 0, 0,  0.2],
                       [3.0, 0, 0, -0.3],
                       [4.5, 0, 0,  0.4],
                       [6.0, 0, 0, -0.5]])

CHAIN_XYZQ_LIST = [[0.0, 0, 0, -0.1],
                   [1.5, 0, 0,  0.2],
                   [3.0, 0, 0, -0.3],
                   [4.5, 0, 0,  0.4],
                   [6.0, 0, 0, -0.5]]


def test_flatten_nested():
    """Tests flattening nested lists."""
    assert list(_flatten([1, [2, [3, 4]], (5,)])) == [1, 2, 3, 4, 5]


def test_flatten_strings():
    """Tests flattening lists of strings."""
    assert list(_flatten(["hello", ["there"]])) == ["hello", "there"]
    assert list(_flatten([["master"], "kenobi"])) == ["master", "kenobi"]


def test_normalized_vector():
    """Tests vector normalization."""
    assert np.allclose(normalized_vector([3.0, 0.0, 4.0]), [0.6, 0.0, 0.8])
    assert np.allclose(normalized_vector([1.5, -2.2, 3.0, -0.2]), [0.37348578, -0.54777915,  0.74697157, -0.0497981])


def test_filter_xyzq_numpy():
    """Tests filtering numpy arrays."""
    assert np.allclose(filter_xyzq(CHAIN_XYZQ, [2, 4]), CHAIN_XYZQ[[1, 3]])
    assert np.allclose(filter_xyzq(CHAIN_XYZQ, [2], charges=False), [[1.5, 0, 0]])
    assert np.allclose(filter_xyzq(CHAIN_XYZQ, [2], coordinates=False), [[0.2]])


def test_filter_xyzq_list():
    """Tests filtering lists."""
    assert np.allclose(filter_xyzq(CHAIN_XYZQ_LIST, [2, 4]), [CHAIN_XYZQ_LIST[1], CHAIN_XYZQ_LIST[3]])
    assert np.allclose(filter_xyzq(CHAIN_XYZQ_LIST, [2], charges=False), [[1.5, 0, 0]])
    assert np.allclose(filter_xyzq(CHAIN_XYZQ_LIST, [2], coordinates=False), [[0.2]])


def test_mask_atoms():
    """Tests masking (setting to zero) atomic info in an array."""
    forces = np.arange(15, dtype=float).reshape(5, 3)
    masked = mask_atoms(forces, [1, 3])
    assert np.allclose(masked[[0, 2]], forces[[0, 2]])
    assert np.allclose(masked[[1, 3, 4]], 0.0)


