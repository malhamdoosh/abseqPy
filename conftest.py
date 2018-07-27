import pytest
import numpy


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    doctest_namespace['np'] = numpy


@pytest.fixture(autouse=True)
def add_raises(doctest_namespace):
    doctest_namespace['raises'] = pytest.raises
