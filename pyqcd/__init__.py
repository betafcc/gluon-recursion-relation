from typing import Literal
from dataclasses import dataclass


@dataclass(frozen=True)
class VectorInfo:
    helicity: Literal[-1, 1]
