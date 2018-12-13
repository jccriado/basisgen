from invariants.weights import Weight
from invariants.representations import Irrep

import abc
import collections
import enum
import functools
import itertools
import re
import operator
import numpy as np


class Series(enum.Enum):
    A = 1
    B = 2
    C = 3
    D = 4
    E = 5
    F = 6
    G = 7


class Algebra(metaclass=abc.ABCMeta):
    def _to_semisimple(self):
        pass

    @property
    def level_vector(self):
        pass

    def height(self, weight):
        return sum(map(operator.mul, self.level_vector, weight))


class SimpleAlgebra(Algebra):
    def __init__(self, series, rank):
        self.series = series
        self.rank = rank

    def __str__(self):
        series_str = {
            Series.A: 'A',
            Series.B: 'B',
            Series.C: 'C',
            Series.D: 'D',
            Series.E: 'E',
            Series.F: 'F',
            Series.G: 'G'
        }[self.series]

        return "{series}{rank}".format(series=series_str, rank=self.rank)

    __repr__ = __str__

    def __hash__(self):
        return hash((self.series, self.rank))

    def __eq__(self, other):
        if not isinstance(other, SimpleAlgebra):
            return False
        else:
            return self.series == other.series and self.rank == other.rank

    def __add__(self, other):
        return self._to_semisimple() + other._to_semisimple()

    def _to_semisimple(self):
        return SemisimpleAlgebra([self])

    @property
    def cartan_matrix(self):
        def default_element(i, j):
            return {0: 2, 1: -1}.get(abs(i - j), 0)

        def generic_cartan_matrix(exceptional_elements):
            return [
                [
                    exceptional_elements.get((i, j), default_element(i, j))
                    for j in range(self.rank)
                ]
                for i in range(self.rank)
            ]

        exceptional_elements = {
            Series.A: {},
            Series.B: {(self.rank - 2, self.rank - 1): -2},
            Series.C: {(self.rank - 1, self.rank - 2): -2},
            Series.D: {
                (self.rank - 3, self.rank - 1): -1,
                (self.rank - 1, self.rank - 3): -1,
                (self.rank - 2, self.rank - 1): 0,
                (self.rank - 1, self.rank - 2): 0
            }
        }[self.series]

        return generic_cartan_matrix(exceptional_elements)

    @property
    def simple_roots(self):
        return list(map(Weight, self.cartan_matrix))

    @property
    @functools.lru_cache(maxsize=None)
    def highest_root(self):
        if self.series == Series.A:
            return Weight(
                [1] + [0] * (self.rank - 2) + [1] if self.rank > 1 else [2]
            )
        elif self.series == Series.B or self.series == Series.D:
            return Weight(
                [0, 1] + [0] * (self.rank - 2)
            )
        elif self.series == Series.C:
            return Weight(
                [2] + [0] * (self.rank - 1)
            )

    @property
    @functools.lru_cache(maxsize=None)
    def level_of_simple_roots(self):
        return {
            Series.A: self.rank - 1,
            Series.B: 2 * self.rank - 2,
            Series.C: 2 * self.rank - 2,
            Series.D: 2 * self.rank - 4
        }[self.series]

    @property
    @functools.lru_cache(maxsize=None)
    def positive_roots(self):
        roots = Irrep(self, self.highest_root).weights_by_level
        return list(
            itertools.chain.from_iterable(
                roots[level]
                for level in range(self.level_of_simple_roots + 1)
            )
        )

    @property
    @functools.lru_cache(maxsize=None)
    def sum_of_positive_roots(self):
        return Weight([1] * self.rank)

    @property
    @functools.lru_cache(maxsize=None)
    def metric(self):
        def build_matrix(element, size):
            return [[element(i, j) for j in range(size)] for i in range(size)]

        n = self.rank

        def element_A(i, j):
            return min(i + 1, j + 1) * (n - max(i, j)) / (n + 1)

        def element_B(i, j):
            if i + 1 == n and j + 1 == n:
                return n/4
            elif i + 1 == n:
                return (j + 1)/2
            elif j + 1 == n:
                return (i + 1)/2
            else:
                return min(i + 1, j + 1)

        def element_C(i, j):
            return min(i + 1, j + 1)/2

        def element_D(i, j):
            if i + 1 > n - 2 and j + 1 > n - 2:
                return n/4 if i == j else (n - 2)/4
            elif i + 1 > n - 2:
                return (j + 1)/2
            elif j + 1 > n - 2:
                return (i + 1)/2
            else:
                return min(i + 1, j + 1)

        element = {
            Series.A: element_A,
            Series.B: element_B,
            Series.C: element_C,
            Series.D: element_D
        }[self.series]

        return build_matrix(element, n)

    @functools.lru_cache(maxsize=None)
    def scalar_product(self, first_weight, second_weight):
        return sum(
            first_component * metric_element * second_component
            for first_component, row in zip(first_weight, self.metric)
            for metric_element, second_component in zip(row, second_weight)
        )

    @functools.lru_cache(maxsize=None)
    def norm_squared(self, weight):
        return self.scalar_product(weight, weight)

    @property
    @functools.lru_cache(maxsize=None)
    def level_vector(self):
        if self.series == Series.A:
            return [(self.rank - i) * (i + 1) for i in range(self.rank)]
        elif self.series == Series.B:
            return (
                [
                    (2*self.rank - i) * (i + 1)
                    for i in range(self.rank - 1)
                ]
                + [int(self.rank * (self.rank + 1)/2)]
            )
        elif self.series == Series.C:
            return [
                (2*self.rank - i - 1) * (i + 1)
                for i in range(self.rank)
            ]
        elif self.series == Series.D:
            return (
                [
                    (2*self.rank - i - 2) * (i + 1)
                    for i in range(self.rank - 2)
                ]
                + [int(self.rank * (self.rank - 1)/2)] * 2
            )


class SemisimpleAlgebra(collections.Iterable, Algebra):
    def __init__(self, simple_algebras):
        self.simple_algebras = simple_algebras
        self.rank = sum(
            simple_algebra.rank for simple_algebra in simple_algebras
        )

    def __str__(self):
        return " + ".join(map(str, self.simple_algebras))

    __repr__ = __str__

    def __hash__(self):
        return hash(tuple(self.simple_algebras))

    def __eq__(self, other):
        if not isinstance(other, SemisimpleAlgebra):
            return False
        else:
            return all(
                self_algebra == other_algebra
                for self_algebra, other_algebra in zip(self, other)
            )

    def __add__(self, other):
        return SemisimpleAlgebra(
            self.simple_algebras + other._to_semisimple().simple_algebras
        )

    def __iter__(self):
        return iter(self.simple_algebras)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return SemisimpleAlgebra(self.simple_algebras[key])
        else:
            return self.simple_algebras[key]

    @staticmethod
    def from_strings(*strings):
        return SemisimpleAlgebra(list(map(SimpleAlgebra.from_string, strings)))

    def _to_semisimple(self):
        return self

    def split_weight(self, weight):
        ranks = [simple_algebra.rank for simple_algebra in self]

        split_weight = []
        for rank in ranks:
            split_weight.append(weight[:rank])
            weight = weight[rank:]

        return split_weight

    def join_weights(self, weights):
        return Weight(list(itertools.chain.from_iterable(weights)))

    @property
    @functools.lru_cache(maxsize=None)
    def level_vector(self):
        return Weight(
            list(
                itertools.chain.from_iterable(
                    simple_algebra.level_vector
                    for simple_algebra in self
                )
            )
        )
