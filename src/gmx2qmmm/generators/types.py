from typing import Tuple, Literal

import numpy as np


PointChargeField = np.ndarray[Tuple[int, Literal[4]], np.float64]
# array of shape (n_atoms, 4) representing x, y, z, charge