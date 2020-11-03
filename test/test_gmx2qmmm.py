import gmx2qmmm
import pytest


@pytest.mark.parametrize("nested,flattened", [([1, ["a", "ab"]], [1, "a", "ab"])])
def test_flatten(nested, flattened):
    assert flattened == list(gmx2qmmm._flatten(nested))
