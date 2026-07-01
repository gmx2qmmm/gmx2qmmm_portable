import pathlib
from abc import ABC, abstractmethod
from collections.abc import Iterable
from typing import List, Optional, Union, TextIO

from gmx2qmmm.generators.types import PointChargeField


class PCFGenerator(ABC):
    """Abstract base class for point charge field generators"""

    @abstractmethod
    def generate(self) -> PointChargeField:
        """Generate a point charge field"""


def dump_field(
    fp: TextIO,
    field: PointChargeField,
    *,
    annotations: Optional[Iterable[Union[str, None]]] = None,
    precision: int = 10,
    legacy: bool = False,
) -> None:
    """Dump a point charge field to a file-like object in the format of a PCF file"""

    if annotations is not None:
        annotations = list(annotations)
        if len(annotations) != field.shape[0]:
            raise ValueError(
                f"Length of annotations {len(annotations)} must match the number of charges in the field {field.shape[0]}"
            )
    else:
        annotations = repeat(None)

    if legacy:
        _dump_field_legacy(fp, field, annotations=annotations, precision=precision)
        return

    charge_str = f"{{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}}"
    for charge, annotation in zip(field, annotations):
        if annotation is not None:
            annotation_str = f"  # {annotation}"
        else:
            annotation_str = ""

        fp.write(charge_str.format(*charge) + annotation_str + "\n")


def _dump_field_legacy(
    fp: TextIO,
    field: PointChargeField,
    annotations: Iterable[Union[str, None]],
    precision: int = 10,
) -> None:

    charge_str = f"{{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}} {{:{precision + 4}.{precision}f}}"
    for charge, annotation in zip(field, annotations):
        if annotation in {"QM", "M1"}:
            fp.write("QM\n")
            continue
        fp.write(charge_str.format(*charge) + "\n")
    fp.write("$end\n")


def load_field(fp: TextIO) -> PointChargeField:
    """Load a point charge field from a file-like object in the format of a PCF file"""

    charges = []
    for line in fp:
        line = line.strip()
        if not line or line.startswith("$end"):
            continue

        if line.startswith("QM"):
            # Legacy handling for "QM" lines
            line = "0.0 0.0 0.0 0.0"

        # Strip off annotation if present (anything after a "#" character)
        try:
            data, annotation = line.split("#", 1)
            annotation = annotation.strip()
        except ValueError:
            data = line
            annotation = None

        parts = data.split()
        if len(parts) != 4:
            raise ValueError(f"Invalid line in PCF file: {line}")
        try:
            parsed = [float(part) for part in parts[:4]]
        except ValueError as e:
            raise ValueError(f"Invalid numeric value in line: {line}") from e

        charges.append(parsed)

    return np.array(charges)


def load_field_legacy(fp: Union[StrPath, TextIO]) -> List[List[Union[float, str]]]:
    """Load a point charge field from a file-like object in the legacy format of a PCF file"""

    if isinstance(fp, (str, pathlib.Path)):
        with open(fp, "r") as f:
            return load_field_legacy(f)

    charges = []
    for line in fp:
        line = line.strip()
        if not line or line.startswith("$end"):
            continue

        if line.startswith("QM"):
            charges.append(["QM"])
            continue

        parts = line.split()
        if len(parts) != 4:
            raise ValueError(f"Invalid line in PCF file: {line}")
        try:
            parsed = [float(part) for part in parts]
        except ValueError as e:
            raise ValueError(f"Invalid numeric value in line: {line}") from e

        charges.append(parsed)

    return charges
