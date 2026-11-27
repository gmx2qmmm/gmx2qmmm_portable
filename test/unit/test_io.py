import pytest
from gmx2qmmm.input_handler import FileReader
from gmx2qmmm.logging_utils import Output


def test_read_file(tmp_path):
    """Test basic file parsing."""
    tmp_param_file = tmp_path / "params.txt"
    tmp_param_file.write_text("# comment\nincludeprio=0\n\njobname=test\ncores=4\n"
                 "stepsize=0.5\nextra=\ninneratomslist=[]\n")
    dictionary = FileReader.read_file(None, str(tmp_param_file))
    assert dictionary["jobname"] == "test"       # checks reading strings
    assert dictionary["cores"] == 4              # checks integer reading and conversion
    assert dictionary["stepsize"] == "0.5"       # checks floats reading (kept as strings)
    assert "extra" not in dictionary             # checks that empty values are dropped
    assert dictionary["inneratomslist"] == []    # checks behavior for empty lists


def test_read_file_error_catch(tmp_path):
    """A missing parameter file forces the program to exit with code 1."""
    with pytest.raises(SystemExit) as exc:
        FileReader.read_file(None, str(tmp_path / "non_existing_file.txt"))
    assert exc.value.code == 1


def test_output_files_format(tmp_path):
    """Tests file appending."""    
    out = Output(tmp_path)
    out.oenergy_append(5, -1.0, 2.5, 0.125, 1.625)
    out.oforces_append(7, [[0.1, -0.2, 0.3]])
    
    assert (tmp_path / "oenergy.txt").read_text().splitlines() == [
        "Step\tQM\t\tMM\t\tLink\t\tTotal",
        "5\t-1.000000\t2.500000\t0.125000\t1.625000",
    ]
    assert (tmp_path / "oforces.txt").read_text() == \
        "Step 7\n1 0.10000000 -0.20000000 0.30000000\n\n"

    # let's verify that it indeed appends things
    out.oforces_append(8, [[0.5, 0.25, -0.125]])

    assert (tmp_path / "oforces.txt").read_text() == (
        "Step 7\n"
        "1 0.10000000 -0.20000000 0.30000000\n"
        "\n"
        "Step 8\n"
        "1 0.50000000 0.25000000 -0.12500000\n"
        "\n"
    )


@pytest.mark.parametrize("prio, expected", [(0, "BP86"), (1, "B3LYP")])
def test_read_file_priority_rules(tmp_path, prio, expected):
    """Tests !INCLUDE option of reading parameters from another file. 
    includeprio=0: main file has a higher priority;
    includeprio=1: included file has a higher priority."""
    lib = tmp_path / "lib.txt"
    lib.write_text("method=B3LYP\nbasis=STO-3G\n")
    param_file = tmp_path / "params.txt"
    param_file.write_text(f"includeprio={prio}\nmethod=BP86\n!INCLUDE={lib}\n")
    d = FileReader.read_file(None, str(param_file))
    assert d["method"] == expected
    assert d["basis"] == "STO-3G"


