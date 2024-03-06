from dataclasses import dataclass
from typing import Literal

import numpy as np

Number = complex
type Helicity = Literal[-1, 1]
type Vector = np.ndarray
type Matrix = np.ndarray


@dataclass
class Gluon:
    helicity: Helicity
    momentum: np.ndarray


# dirac matrices
Gamma = (
    np.array(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, -1, 0],
            [0, 0, 0, -1],
        ],
        dtype=Number,
    ),
    np.array(
        [
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [0, -1, 0, 0],
            [-1, 0, 0, 0],
        ],
        dtype=Number,
    ),
    np.array(
        [
            [0, 0, 0, -1j],
            [0, 0, 1j, 0],
            [0, 1j, 0, 0],
            [-1j, 0, 0, 0],
        ],
        dtype=Number,
    ),
    np.array(
        [
            [0, 0, 1, 0],
            [0, 0, 0, -1],
            [-1, 0, 0, 0],
            [0, 1, 0, 0],
        ],
        dtype=Number,
    ),
)


def mp(a: Vector, b: Vector) -> Number:
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def mp2(a: Vector):
    return mp(a, a)


def braket(q: Vector, k: Vector):
    return np.conj(q) @ k


def box(q: Vector, k: Vector):
    # return braket(q, k) - braket(k, q)
    return q @ np.conj(k)


def fish(q: Vector, g: Matrix, k: Vector):
    return q @ (g @ k)


def backfish(q: Vector, g: Matrix, k: Vector):
    return np.conj(q) @ (g @ k)


@dataclass
class Process:
    gluons: list[Gluon]

    def current(self, gis: list[int], xi: int | None = None):
        if xi is None:
            if len(gis) == 1:
                return self.polarization(gis[0])
            elif len(gis) == 2:
                return np.array(
                    [
                        self.current(gis, 0),
                        self.current(gis, 1),
                        self.current(gis, 2),
                        self.current(gis, 3),
                    ],
                    dtype=Number,
                )
            elif len(gis) == 3:
                return (1.0 / mp2(self.kappa(gis))) * (
                    self.sb([0], [1, 2]) + self.sb([0, 1], [2]) + self.cb([0], [1], [2])
                )
            else:
                factor = 1.0 / mp2(self.kappa(gis))
                n = len(gis)

                sb_acc = np.zeros(4, dtype=Number)
                for m in range(0, n - 1):
                    sb_acc += self.sb(gis[: m + 1], gis[m + 1 :])

                cb_acc = np.zeros(4, dtype=Number)
                for m in range(0, n - 2):
                    for k in range(m + 1, n - 1):
                        cb_acc += self.cb(
                            gis[: m + 1], gis[m + 1 : k + 1], gis[k + 1 :]
                        )

                return factor * (sb_acc + cb_acc)
        else:
            if len(gis) == 1:
                return self.polarization(gis[0], xi)
            elif len(gis) == 2:
                return 1.0 / mp2(self.kappa(gis)) * self.sb([0], [1], xi)
            elif len(gis) == 3:
                return (1.0 / mp2(self.kappa(gis))) * (
                    self.sb([0], [1, 2], xi) + self.sb([0, 1], [2], xi)
                )
            else:
                raise NotImplementedError

    def kappa(self, gis: list[int]):
        return sum([self.gluons[i].momentum for i in gis])

    def sb(self, xs: list[int], ys: list[int], xi: int | None = None):
        if xi is None:
            return np.array(
                [
                    self.sb(xs, ys, 0),
                    self.sb(xs, ys, 1),
                    self.sb(xs, ys, 2),
                    self.sb(xs, ys, 3),
                ],
                dtype=Number,
            )
        else:
            return (
                (2 * self.kappa(ys)) @ (self.current(xs) * self.current(ys, xi))
                - (2 * self.kappa(xs)) @ (self.current(ys) * self.current(xs, xi))
                + (self.kappa(xs) - self.kappa(ys))[xi]
                * (self.current(xs) @ self.current(ys))
            )

    def cb(self, a: list[int], b: list[int], c: list[int], xi: int | None = None):
        if xi is None:
            return np.array(
                [
                    self.cb(a, b, c, 0),
                    self.cb(a, b, c, 1),
                    self.cb(a, b, c, 2),
                    self.cb(a, b, c, 3),
                ],
                dtype=Number,
            )
        else:
            return self.current(a) @ (
                self.current(c) * self.current(b, xi)
                - self.current(b) * self.current(c, xi)
            ) - self.current(c) @ (
                self.current(b) * self.current(a, xi)
                - self.current(a) * self.current(b, xi)
            )

    def polarization(self, gi: int, xi: int | None = None):
        if xi is None:
            return np.array(
                [
                    self.polarization(gi, 0),
                    self.polarization(gi, 1),
                    self.polarization(gi, 2),
                    self.polarization(gi, 3),
                ],
                dtype=Number,
            )
        else:
            gm = Gamma[xi]
            q = self.auxiliar(gi)
            k = self.gluons[gi].momentum

            if self.gluons[gi].helicity == 1:
                return (-1 / (2**0.5)) * (backfish(q, gm, k) / braket(q, k))
            elif self.gluons[gi].helicity == -1:
                return (1 / (2**0.5)) * (fish(q, gm, k) / box(q, k))
            else:
                raise ValueError

    def auxiliar(self, gi: int):
        return self.gluons[(gi + 1) % len(self.gluons)].momentum
