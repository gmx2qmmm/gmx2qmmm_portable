"""Regression test for PCF generation"""

import shutil
import pathlib
from dataclasses import dataclass
from typing import List

import numpy as np
import pytest

from gmx2qmmm.generators.pcf import GeneratePCF, PCFGeneratorShift, dump_field


@dataclass
class MockSystem:
    array_xyzq_current: np.ndarray
    list_atoms_qm: List[int]
    list_atoms_m1: List[int]
    list_atoms_m2: List[List[int]]


@pytest.fixture
def system(pcf_input_data):
    return MockSystem(
        array_xyzq_current=pcf_input_data["xyzq"],
        list_atoms_qm=pcf_input_data["qm_atoms"],
        list_atoms_m1=pcf_input_data["m1_atoms"],
        list_atoms_m2=pcf_input_data["m2_atoms"],
    )


@pytest.fixture
def input_dict(pcf_input_data):
    return {
        "jobname": pcf_input_data["jobname"],
        "charge": pcf_input_data["charge"],
    }


def test_pcf_old(system, input_dict):
    """Test that PCF generation produces expected output"""

    module_path = pathlib.Path(__file__)
    work_dir = module_path.parent / module_path.stem
    work_dir.mkdir(exist_ok=True)

    pcf = GeneratePCF(
        input_dict=input_dict,
        system=system,
        topology=None,  # Not used
        work_dir=work_dir,
    )

    output_file = work_dir / "test.pointcharges"
    output_file = output_file.rename(work_dir / "test.pointcharges.current")
    reference_output_file = work_dir / "test.pointcharges.ref"

    if not reference_output_file.exists():
        shutil.move(output_file, reference_output_file)
        assert False, (
            "Reference output file did not exist, created from current output. Please verify that the generated PCF is correct and commit the reference output file."
        )
    else:
        with open(output_file) as f_out, open(reference_output_file) as f_ref:
            output_lines = f_out.readlines()
            reference_lines = f_ref.readlines()
            assert output_lines == reference_lines, (
                "Generated PCF does not match reference output"
            )

    print(
        f"delta={np.linalg.norm(pcf.new_sum - pcf.target_sum)} after {pcf.iterations_} iterations"
    )


def test_pcf_new(system, input_dict):
    """Test that PCF generation produces expected output"""

    module_path = pathlib.Path(__file__)
    work_dir = module_path.parent / module_path.stem
    work_dir.mkdir(exist_ok=True)

    pcf = PCFGeneratorShift.from_system(system, charge=input_dict["charge"])
    field = pcf.generate()

    output_file = work_dir / "test.pointcharges.new"
    reference_output_file = work_dir / "test.pointcharges.ref"

    with open(output_file, "w") as fp:
        dump_field(fp, field, annotations=pcf.annotations)

    if not reference_output_file.exists():
        assert False, (
            "Reference output file does not exist. Please run test_pcf_current first to generate it."
        )
    else:
        with open(output_file) as fp, open(reference_output_file) as fp_ref:
            output_lines = fp.readlines()
            reference_lines = fp_ref.readlines()

        correction_lines = []
        correction_lines_ref = []
        for i, (line_out, line_ref) in enumerate(zip(output_lines, reference_lines), 1):
            if line_ref.strip().startswith("QM"):
                continue

            if line_ref.strip() == "$end":
                break

            line_out, annotation = line_out.split("#", 1)

            x, y, z, q = map(float, line_out.split())
            x_ref, y_ref, z_ref, q_ref = map(float, line_ref.split())

            if annotation.strip() == "CORRECTION":
                correction_lines.append(np.array([x, y, z, q]))
                correction_lines_ref.append(np.array([x_ref, y_ref, z_ref, q_ref]))
                continue

            assert np.allclose([x, y, z, q], [x_ref, y_ref, z_ref, q_ref], atol=1e-6), (
                f"Generated PCF line {i} does not match reference output:\n"
                f"Output: {line_out}\n"
                f"Reference: {line_ref}"
            )

        for a in correction_lines:
            assert any(np.allclose(a, b, atol=1e-3) for b in correction_lines_ref), (
                f"Generated correction charge {a} does not match any reference correction charge"
            )

        print(
            f"delta={np.linalg.norm(pcf.current_field_vector_ - pcf.target_field_vector_)} after {pcf.iterations_} iterations"
        )
