import numpy as np

from gmx2qmmm.generators.geometry import make_xyzq


def test_make_xyzq():

    coords  = [0.123, 0.456, 0.789, 0.321, 0.654, 0.987]
    charges = [-0.57, +1.57]
    ref_xyzq = np.array([
        [0.123, 0.456, 0.789, -0.57],
        [0.321, 0.654, 0.987, +1.57],
        ])

    res_xyzq = make_xyzq(coords,charges)

    assert np.allclose(res_xyzq, ref_xyzq)
