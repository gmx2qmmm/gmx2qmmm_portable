import pytest
import numpy as np
from pathlib import Path

from gmx2qmmm.app import App

# helper function
def clean_up():
    current_path = Path(__file__).resolve().parent
    work_dir = current_path / 'test_files'

    files_to_delete = [ 'test.pointcharges'
                      , 'logfile.log'
                      , 'test.qmmm.top'
                      , 'test.boxlarge.g96'
                      , 'test.gjf'
                      , 'test.mdp'
                      , 'test.gjf.log'
                      , 'test.qmmm.top.ndx'
                      , 'test.gmx.log']

    for i in files_to_delete:
        file = work_dir / i
        file.unlink(missing_ok=True)

   
def test_app():
    test_path = Path(__file__).resolve().parent / 'test_files'
    ref_path = Path(__file__).resolve().parent / 'ref_files'
    
    files_to_check = ['test.pointcharges', 'test.qmmm.top', 'test.qmmm.top.ndx'] 
    app = App(parameters=str(test_path / 'params.txt'), work_dir=str(test_path)) 
   
    # check whether the path to logfile.log is correct
    assert(app.logfile == test_path / 'logfile.log')
     
    # check line by line whether created files contain correct content 
    for file_name in files_to_check: 
        with open(test_path / file_name, 'r') as a, open(ref_path / file_name, 'r') as b:
            for line_a, line_b in zip(a,b):
                assert line_a.rstrip() == line_b.rstrip()
    clean_up()


if __name__ == '__main__':
    pytest.main([__file__]) 
