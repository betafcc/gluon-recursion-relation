from random import randint

from tests.util import random_process
from numpy.testing import assert_array_almost_equal


def test_random_process():
    for _ in range(10):
        process = random_process(randint(3, 10))
        momenta = [g.momentum for g in process.gluons]

        assert_array_almost_equal(sum(momenta[:2]), sum(momenta[2:]), decimal=10)
