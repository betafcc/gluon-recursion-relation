import numpy as np
from numpy.testing import assert_array_almost_equal

from pyqcd.process import Gluon, Process
from tests.util import random_process


process = random_process(5)


def test_kappa_dot_current():
    assert (process.kappa([0, 1]) @ process.current([0, 1])).all(0)


def test_210_012():
    assert_array_almost_equal(process.current([2, 1, 0]), process.current([0, 1, 2]))


def test_perm_0():
    assert_array_almost_equal(
        process.current([0, 1, 2])
        + process.current([1, 2, 0])
        + process.current([2, 0, 1]),
        np.array([0, 0, 0, 0], dtype=complex),
    )
