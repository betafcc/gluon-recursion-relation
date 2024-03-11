import math
import random
from typing import Sequence

import numpy as np


def random_angle():
    return random.uniform(0, 2 * math.pi)


def random_momentum(mass: float, energy: float):
    theta = random_angle()
    cos_theta = math.cos(theta)
    sin_theta = (1 - cos_theta * cos_theta) ** 0.5

    phi = 2 * math.pi * random.random()

    pmod = math.sqrt(energy * energy - mass * mass)
    return np.array(
        [
            pmod * sin_theta * math.cos(phi),
            pmod * sin_theta * math.sin(phi),
            pmod * cos_theta,
            energy,
        ]
    )


# def invariant_mass(cme: float, masses: Sequence[float]):
#     n = len(masses)
#     invariant_m = [0.0] * n
#     mass_total = sum(masses)
#     random_numbers = sorted([random.random() for _ in range(n - 1)])

#     acc = 0.0
#     for i in range(len(masses) - 1):
#         acc += masses[i]
#         invariant_m[i] = acc + random_numbers[i] * (cme**0.5 - mass_total)

#     invariant_m[-1] = cme**0.5

#     return invariant_m


def invariant_mass(cme: float, masses: Sequence[float]):
    from random import random
    from itertools import accumulate

    diff = cme**0.5 - sum(masses)
    return [
        m + seed * diff
        for m, seed in zip(
            accumulate(masses),
            sorted(random() for _ in range(len(masses) - 1)),
        )
    ] + [cme**0.5]


def lorentz_transform(v: np.ndarray, boost: np.ndarray) -> np.ndarray:
    beta2 = boost[0] ** 2 + boost[1] ** 2 + boost[2] ** 2
    gamma = 1.0 / (1 - beta2) ** 0.5
    factor1 = (gamma - 1.0) / beta2
    factor2 = v[0] * boost[0] + v[1] * boost[1] + v[2] * boost[2]

    return np.array(
        [
            v[0] + boost[0] * factor1 * factor2 + gamma * boost[0] * v[3],
            v[1] + boost[1] * factor1 * factor2 + gamma * boost[1] * v[3],
            v[2] + boost[2] * factor1 * factor2 + gamma * boost[2] * v[3],
            gamma * (v[3] + factor2),
        ]
    )


def initial_pair(cme: float):
    z = (cme**0.5) / 2

    return (
        np.array([0, 0, z, z]),
        np.array([0, 0, -z, z]),
    )


def phase_space(cme: float, masses: Sequence[float]) -> Sequence[np.ndarray]:
    pa, pb = initial_pair(cme)

    invariant_mass_values = invariant_mass(cme, masses)

    E2: float = (
        invariant_mass_values[1] * invariant_mass_values[1]
        + masses[1] * masses[1]
        - masses[0] * masses[0]
    ) / (2 * invariant_mass_values[1])

    E1: float = invariant_mass_values[1] - E2

    p2 = random_momentum(masses[1], E2)
    p1 = np.array([-p2[0], -p2[1], -p2[2], E1])

    final_state_particles = [p1, p2]

    for i in range(3, len(masses) + 1):
        Ek = (
            invariant_mass_values[i - 1] * invariant_mass_values[i - 1]
            + invariant_mass_values[i - 2] * invariant_mass_values[i - 2]
            - masses[i - 1] * masses[i - 1]
        ) / (2 * invariant_mass_values[i - 1])

        E = invariant_mass_values[i - 1] - Ek

        p = random_momentum(masses[i - 1], E)
        final_state_particles.append(p)

        boost = np.array([-p[0] / Ek, -p[1] / Ek, -p[2] / Ek])
        for j in range(i - 1):
            final_state_particles[j] = lorentz_transform(
                final_state_particles[j], boost
            )

    return [pa, pb, *final_state_particles]


def random_process(n: int, cme: float | None = None):
    from pyqcd.process import Process, Gluon, Number

    cme = cme if cme is not None else random.uniform(200, 13000)

    if n == 0:
        raise ValueError("n must be greater than 0")

    return Process(
        [
            Gluon(random.choice([-1, 1]), np.array(list(p), dtype=Number))
            for p in phase_space(cme, [0.0] * n)
        ]
    )
