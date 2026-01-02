import pytest
from pathlib import Path
from unittest.mock import patch

from gmx2qmmm.app import App
from gmx2qmmm.jobs import qm, mm

def test_singlepoint():

    current_path = Path(__file__).resolve().parent
    params_file = current_path / "params.txt"

    # Create App instance
    app = App(parameters=str(params_file), work_dir=str(current_path ))

    # Mock the main components
    # mock_singlepoint = mocker.patch.object("gmx2qmmm.app.Singlepoint")
    # mock_qm = mocker.patch.object(mock_singlepoint.class_qm_job, 'run_qm_job')
    # mock_gmx = mocker.patch.object(mock_singlepoint.class_mm_job, 'run_gmx')

    with patch.object(qm.QM_gaussian, 'run_qm_job'), \
         patch.object(mm.MM, 'run_gmx'):

        app.run()



if __name__ == '__main__':
    pytest.main([__file__]) 
