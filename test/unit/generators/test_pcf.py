"""Regression test for PCF generation"""

import shutil
import pathlib
from dataclasses import dataclass
from typing import List

from gmx2qmmm.generators.pcf import GeneratePCF


@dataclass
class MockSystem:
    array_xyzq_current: np.ndarray
    list_atoms_qm: List[int]
    list_atoms_m1: List[int]
    list_atoms_m2: List[List[int]]


def test_pcf_current(pcf_input_data):
    """Test that PCF generation produces expected output for current input data"""

    module_path = pathlib.Path(__file__)
    work_dir = module_path.parent / module_path.stem
    work_dir.mkdir(exist_ok=True)

    system = MockSystem(
        array_xyzq_current=pcf_input_data['xyzq'],
        list_atoms_qm=pcf_input_data['qm_atoms'],
        list_atoms_m1=pcf_input_data['m1_atoms'],
        list_atoms_m2=pcf_input_data['m2_atoms'],
    )

    input_dict = {
        "jobname": pcf_input_data['jobname'],
        "charge": pcf_input_data['charge'],
    }

    pcf = GeneratePCF(
        input_dict=input_dict,
        system=system,
        topology=None,  # Not used
        work_dir=work_dir,
    )

    output_file = work_dir / 'test.pointcharges'
    reference_output_file = work_dir / 'test.pointcharges.ref'

    if not reference_output_file.exists():
        shutil.move(output_file, reference_output_file)
    else:
        with open(output_file) as f_out, open(reference_output_file) as f_ref:
            output_lines = f_out.readlines()
            reference_lines = f_ref.readlines()
            assert output_lines == reference_lines, "Generated PCF does not match reference output"
