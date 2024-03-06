from dataclasses import dataclass
from typing import Literal

import numpy as np

Number = complex
type Helicity = Literal[-1, 1]
type Vector = np.ndarray
type Matrix = np.ndarray


@dataclass
class Gluon:
    polarization: np.ndarray
    momentum: np.ndarray


def mp(a: Vector, b: Vector) -> Number:
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def mp2(a: Vector):
    return mp(a, a)


@dataclass
class Process:
    gluons: list[Gluon]

    def current(self, gis: list[int], xi: int | None = None):
        if xi is None:
            if len(gis) == 1:
                return self.gluons[gis[0]].polarization
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
                return self.gluons[gis[0]].polarization[xi]
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
