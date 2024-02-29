from __future__ import annotations

from typing import Literal, Union
from dataclasses import dataclass


Number = Union[float, int, complex]


class Vector:
    def __init__(self, x0: Number, x1: Number, x2: Number, x3: Number):
        self.vector = x0, x1, x2, x3

    def __repr__(self):
        return "⟨ " + " ".join(map(str, self.vector)) + " ⟩"

    def __add__(self, other: Vector) -> Vector:
        return Vector(*[x + y for x, y in zip(self.vector, other.vector)])

    def __sub__(self, other: Vector) -> Vector:
        return Vector(*[x - y for x, y in zip(self.vector, other.vector)])

    def __mul__(self, other: Number) -> Vector:
        return Vector(*[x * other for x in self.vector])

    def __rmul__(self, other: Number) -> Vector:
        return self * other

    def __truediv__(self, other: Number) -> Vector:
        return Vector(*[x / other for x in self.vector])

    def __neg__(self) -> Vector:
        return Vector(*[-x for x in self.vector])

    def __eq__(self, other: Vector) -> bool:
        return self.vector == other.vector

    def __getitem__(self, index: Literal[0, 1, 2, 3]) -> Number:
        return self.vector[index]

    def __iter__(self):
        return iter(self.vector)

    def __len__(self):
        return len(self.vector)

    def __matmul__(self, other: Vector) -> Number:
        return sum(x * y for x, y in zip(self.vector, other.vector))

    def __abs__(self) -> Number:
        return sum(x * x for x in self.vector) ** 0.5

    # Minkoswi inner product
    def minkowski(self, other: Vector) -> Number:
        return (
            self[0] * other[0]
            - self[1] * other[1]
            - self[2] * other[2]
            - self[3] * other[3]
        )

    # Lorentz transformation
    def boost(self, beta: Vector) -> Vector:
        gamma = (1 - abs(beta) ** 2) ** -0.5
        return Vector(
            gamma * (self[0] - beta @ self[1:]),
            -gamma * (beta @ self[1:])
            + self[0] * beta[0]
            + self[1] * beta[1]
            + self[2] * beta[2],
            -gamma * (beta @ self[1:])
            + self[0] * beta[0]
            + self[1] * beta[1]
            + self[2] * beta[2],
            -gamma * (beta @ self[1:])
            + self[0] * beta[0]
            + self[1] * beta[1]
            + self[2] * beta[2],
        )


@dataclass
class Gluon:
    helicity: Literal["+", "-"]
    momentum: Vector


@dataclass
class Process:
    gluons: list[Gluon]

    # current recursive calculation for n-gluons
    def amplitude(self) -> complex:
        if len(self.gluons) == 2:
            return self.gluons[0].momentum.minkowski(self.gluons[1].momentum)
        else:
            return sum(
                self.gluons[0].momentum.minkowski(gluon.momentum)
                * Process(gluons=self.gluons[1:]).amplitude()
                for gluon in self.gluons[1:]
            )
