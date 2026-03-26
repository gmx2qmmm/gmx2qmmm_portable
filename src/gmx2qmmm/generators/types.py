from typing import Tuple, Literal

import numpy as np


Vector = np.ndarray[Tuple[int], np.float64]
Vectors = np.ndarray[Tuple[int, int], np.float64]
# General-purpose array types for vectors and collections of vectors

Vector3D = np.ndarray[Tuple[Literal[3]], np.float64]
Vectors3D = np.ndarray[Tuple[int, Literal[3]], np.float64]
# Vector(s) in 3D spaces

Point = Vector3D
Points = Vectors3D
# Technically vector(s) but specifically for particle coordinates in 3D Cartesian space (x, y, z)

PointCharge = np.ndarray[Tuple[Literal[4]], np.float64]
PointChargeField = np.ndarray[Tuple[int, Literal[4]], np.float64]
# Array of shape (n_atoms, 4) representing x, y, z, charge