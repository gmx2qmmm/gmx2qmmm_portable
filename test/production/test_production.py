import pytest
import numpy as np
from pathlib import Path

from gmx2qmmm.app import App
from gmx2qmmm.jobs import qm, mm

def test_singlepoint(mocker):

    current_path = Path(__file__).resolve().parent
    work_dir = current_path / 'test_files'
    params_file = work_dir / 'params.txt'

    # Create App instance
    app = App(parameters=str(params_file), work_dir=str(work_dir))

    mocker.patch.object(qm.QM_gaussian, 'run_qm_job')
    mocker.patch.object(mm.MM, 'run_gmx')
    mocker.patch.object(mm.MM, 'call_mm_forces')
    mocker.patch.object(mm.MM, 'call_mm_energy')
    mocker.patch.object(mm.MM, 'make_gmx_inp')

    # Run app
    app.run()

    # compare energies and forces
    oe_file = work_dir / 'oenergy.txt'
    energies = [float(i) for i in open(oe_file).readlines()[-1].split()]
    ref_energies = [float(i) for i in open('test/production/ref_output/oenergy.txt').readlines()[-1].split()]
    assert np.all(np.isclose(energies, ref_energies))
    
    of_file = work_dir / 'oforces.txt'
    forces = [float(u) for i in [o.split() for o in open(of_file).readlines()[1:-1]] for u in i]
    ref_forces = [float(i) for line in open('test/production/ref_output/oforces.txt').readlines()[1:] for i in line.split()]
    assert np.all(np.isclose(forces, ref_forces))
    
if __name__ == '__main__':
    pytest.main([__file__]) 
