import pytest
import numpy as np
from pathlib import Path

from gmx2qmmm.app import App
from gmx2qmmm.jobs import qm, mm

def test_singlepoint(mocker):

    current_path = Path(__file__).resolve().parent
    work_dir = current_path / 'test_files'
    ref_path = Path(__file__).resolve().parent / 'ref_output'
    
    # Create App instance
    app = App(parameters=str(work_dir / 'params.txt'), work_dir=str(work_dir))

    # check whether the path to logfile.log is correct
    assert(app.logfile == work_dir / 'logfile.log')
 
    # check line by line whether created files contain correct content 
    files_to_check = ['test.pointcharges', 'test.qmmm.top', 'test.qmmm.top.ndx'] 
    for file_name in files_to_check: 
        with open(work_dir / file_name, 'r') as a, open(ref_path / file_name, 'r') as b:
            for line_a, line_b in zip(a,b):
                assert line_a.rstrip() == line_b.rstrip()

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
    ref_energies = [float(i) for i in open(ref_path / 'oenergy.txt').readlines()[-1].split()]
    assert np.all(np.isclose(energies, ref_energies))
    
    of_file = work_dir / 'oforces.txt'
    forces = [float(u) for i in [o.split() for o in open(of_file).readlines()[1:-1]] for u in i]
    ref_forces = [float(i) for line in open(ref_path / 'oforces.txt').readlines()[1:] for i in line.split()]
    assert np.all(np.isclose(forces, ref_forces))

    oe_file.unlink()
    of_file.unlink()
    files_to_delete = [ 'test.pointcharges'
                      , 'logfile.log'
                      , 'test.qmmm.top'
                      , 'test.boxlarge.g96'
                      , 'test.gjf'
                      , 'test.mdp'
                      , 'test.qmmm.top.ndx']

    for i in files_to_delete:
        file = work_dir / i
        file.unlink(missing_ok=True)
    
if __name__ == '__main__':
    pytest.main([__file__]) 
