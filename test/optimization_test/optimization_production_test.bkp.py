import pytest
import numpy as np
from pathlib import Path

import gmx2qmmm.app
import gmx2qmmm.generators.energies
import gmx2qmmm.jobs


def test_optimization(mocker):

    current_path = Path(__file__).resolve().parent
    work_dir = current_path / 'test_files'
    ref_path = Path(__file__).resolve().parent / 'ref_output'
    
    # Create App instance
    app = gmx2qmmm.app.App(parameters=str(work_dir / 'params.txt'), work_dir=str(work_dir))

    # check whether the path to logfile.log is correct
    assert(app.logfile == work_dir / 'logfile.log')
 
    # check line by line whether created files contain correct content 
    files_to_check = ['test.pointcharges', 'test.qmmm.top', 'test.qmmm.top.ndx'] 
    for file_name in files_to_check: 
        with open(work_dir / file_name, 'r') as a, open(ref_path / file_name, 'r') as b:
            for line_a, line_b in zip(a,b):
                assert line_a.rstrip() == line_b.rstrip()

    ref_energies_by_step = {}
    for line in open(ref_path / 'oenergy.txt').readlines()[1:]:
        step, _, _, _, total = line.split()
        ref_energies_by_step[int(step)] = float(total)

    ref_forces_by_step = {}
    current_step = None
    current_forces = []
    for raw_line in open(ref_path / 'oforces.txt').readlines():
        line = raw_line.strip()
        if not line:
            if current_step is not None:
                ref_forces_by_step[current_step] = np.array(current_forces)
                current_step = None
                current_forces = []
            continue
        if line.startswith('Step '):
            current_step = int(line.split()[1])
            current_forces = []
            continue
        _, fx, fy, fz = line.split()
        current_forces.append([float(fx), float(fy), float(fz)])

    if current_step is not None:
        ref_forces_by_step[current_step] = np.array(current_forces)

    def mock_linkcorrenergy(self):
        step = int(self.system.int_step_current)
        qmenergy = float(self.class_qm_job.qmenergy)
        mmenergy = float(self.class_mm_job.mmenergy)
        return qmenergy + mmenergy - ref_energies_by_step[step]

    def mock_linkcorrforces(self):
        step = int(self.system.int_step_current)
        ref_force = ref_forces_by_step[step]
        return np.array(self.class_qm_job.qmforces) + np.array(self.class_mm_job.mmforces) - ref_force

    mocker.patch.object(gmx2qmmm.jobs.qm.QM_gaussian, 'run_qm_job')
    mocker.patch.object(gmx2qmmm.jobs.mm.MM, 'run_gmx')
    mocker.patch.object(gmx2qmmm.jobs.mm.MM, 'call_mm_forces')
    mocker.patch.object(gmx2qmmm.jobs.mm.MM, 'call_mm_energy')
    mocker.patch.object(gmx2qmmm.jobs.mm.MM, 'make_gmx_inp')
    mocker.patch.object(gmx2qmmm.generators.energies.GeneratorEnergies, 'get_linkcorrenergy', autospec=True, side_effect=mock_linkcorrenergy)
    mocker.patch.object(gmx2qmmm.generators.energies.GeneratorForces, 'get_linkcorrforces', autospec=True, side_effect=mock_linkcorrforces)

    # Run app
    app.run()

    # compare energies and forces
    oe_file = work_dir / 'oenergy.txt'
    energies = [float(i) for i in open(oe_file).readlines()[-1].split()]
    ref_energies = [float(i) for i in open(ref_path / 'oenergy.txt').readlines()[-1].split()]
    assert np.all(np.isclose(energies, ref_energies))
    
    of_file = work_dir / 'oforces.txt'

    def parse_force_values(path):
        values = []
        for line in open(path).readlines():
            stripped = line.strip()
            if not stripped or stripped.startswith('Step'):
                continue
            _, fx, fy, fz = stripped.split()
            values.extend([float(fx), float(fy), float(fz)])
        return values

    forces = parse_force_values(of_file)
    ref_forces = parse_force_values(ref_path / 'oforces.txt')
    assert np.all(np.isclose(forces, ref_forces))

    # oe_file.unlink()
    # of_file.unlink()
    files_to_delete = [ 'test.pointcharges', 'test.1.pointcharges', 'test.2.pointcharges', 'test.3.pointcharges', 
                        'test.gjf', 'test.1.gjf', 'test.2.gjf', 'test.3.gjf', 
                        'logfile.log',
                        'test.qmmm.top',
                        'test.boxlarge.g96',
                        'test.mdp',
                        'test.qmmm.top.ndx']

    for i in files_to_delete:
        file = work_dir / i
        file.unlink(missing_ok=True)
    
if __name__ == '__main__':
    pytest.main([__file__])
