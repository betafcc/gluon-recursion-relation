import numpy as np
from numpy.testing import assert_array_almost_equal

from pyqcd.process2 import Gluon, Process


def generate_random_gluon():
    p = np.random.normal(size=3)  # Random spatial components
    p_magnitude = np.linalg.norm(p)
    momentum = np.array([p_magnitude, *p], dtype=complex)  # E=p for massless gluon

    polarization = np.random.normal(size=3)

    polarization -= polarization.dot(p) / p_magnitude**2 * p

    polarization = np.insert(polarization, 0, 0)

    return Gluon(polarization=polarization, momentum=momentum)


gluons = [generate_random_gluon() for _ in range(4)]

process = Process(gluons)


def test_process():

    assert (process.kappa([0, 1]) @ process.current([0, 1])).all(0)

    assert_array_almost_equal(process.current([2, 1, 0]), process.current([0, 1, 2]))

    assert_array_almost_equal(
        process.current([0, 1, 2])
        + process.current([1, 2, 0])
        + process.current([2, 0, 1]),
        np.array([0, 0, 0, 0], dtype=complex),
    )
