import numpy as np
from numpy.testing import assert_array_almost_equal

from pyqcd.process import Gluon, Process


gluons = [
    Gluon(1, np.array([7, 0.2503853218528432, 0.0516972277498482, -6.995329483823019])),
    Gluon(
        -1,
        np.array([17, 1.2068332374603898, 15.309925694446564, -7.290386050648196]),
    ),
    Gluon(
        1,
        np.array([11, -4.839116417149157, -5.099873173941427, -8.460156376272847]),
    ),
    Gluon(1, np.array([18, 5.69161547295297, 6.240582197700386, -15.895302675375119])),
]

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
