from basisgen.weights import Weight

import abc
import collections
import enum
import functools
import itertools
import operator


class Series(enum.Enum):
    A = 1
    B = 2
    C = 3
    D = 4
    E = 5
    F = 6
    G = 7

    def __str__(self):
        return {
            Series.A: 'A',
            Series.B: 'B',
            Series.C: 'C',
            Series.D: 'D',
            Series.E: 'E',
            Series.F: 'F',
            Series.G: 'G'
        }[self]


class Algebra(metaclass=abc.ABCMeta):
    def _to_semisimple(self):
        pass

    @property
    def level_vector(self):
        pass

    def height(self, weight):
        return sum(map(operator.mul, self.level_vector, weight))


class SimpleAlgebra(Algebra):
    class IncorrectRank(Exception):
        def __init__(self, series, rank, rank_bounds):
            series = series
            template_msg = (
                "Unexpected rank {rank} for algebra series {series}. "
                "The rank n must satisfy: {rank_bounds}"
            )
            super().__init__(
                template_msg.format(
                    rank=rank,
                    series=series,
                    rank_bounds=rank_bounds
                )
            )

    def __init__(self, series, rank):
        SimpleAlgebra._check_rank_bounds(series, rank)

        self.series = series
        self.rank = rank

    @staticmethod
    def _check_rank_bounds(series, rank):
        if series == Series.A and rank < 1:
            raise SimpleAlgebra.IncorrectRank(series, rank, "n >= 1")
        if series == Series.B and rank < 2:
            raise SimpleAlgebra.IncorrectRank(series, rank, "n >= 2")
        if series == Series.C and rank < 2:
            raise SimpleAlgebra.IncorrectRank(series, rank, "n >= 2")
        if series == Series.D and rank < 4:
            raise SimpleAlgebra.IncorrectRank(series, rank, "n >= 4")
        if series == Series.E and (rank < 6 or rank > 8):
            raise SimpleAlgebra.IncorrectRank(series, rank, "6 <= n <= 8")
        if series == Series.F and rank != 4:
            raise SimpleAlgebra.IncorrectRank(series, rank, "n = 4")
        if series == Series.G and rank != 2:
            raise SimpleAlgebra.IncorrectRank(series, rank, "n = 2")

    def __str__(self):
        return "{series}{rank}".format(series=self.series, rank=self.rank)

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
            },
            Series.E: {
                (2, self.rank - 1): -1,
                (self.rank - 1, 2): -1,
                (self.rank - 1, self.rank - 2): 0,
                (self.rank - 2, self.rank - 1): 0
            },
            Series.F: {(1, 2): -2},
            Series.G: {(0, 1): -3}
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

        if self.series == Series.B:
            return Weight([0, 1] + [0] * (self.rank - 2))

        if self.series == Series.C:
            return Weight([2] + [0] * (self.rank - 1))

        if self.series == Series.D:
            return Weight([0, 1] + [0] * (self.rank - 2))

        if self.series == Series.E:
            if self.rank == 6:
                return Weight([0, 0, 0, 0, 0, 1])
            elif self.rank == 7:
                return Weight([1, 0, 0, 0, 0, 0, 0])
            elif self.rank == 8:
                return Weight([0, 0, 0, 0, 0, 0, 1, 0])

        if self.series == Series.F:
            return Weight([1, 0, 0, 0])

        if self.series == Series.G:
            return Weight([1, 0])

    @property
    @functools.lru_cache(maxsize=None)
    def level_of_simple_roots(self):
        return {
            Series.A: self.rank - 1,
            Series.B: 2 * self.rank - 2,
            Series.C: 2 * self.rank - 2,
            Series.D: 2 * self.rank - 4,
            Series.E: {6: 10, 7: 16, 8: 28}.get(self.rank, None),
            Series.F: 10,
            Series.G: 4
        }[self.series]

    @property
    @functools.lru_cache(maxsize=None)
    def sum_of_positive_roots(self):
        return Weight([1] * self.rank)

    @property
    @functools.lru_cache(maxsize=None)
    def metric(self):
        if self.series == Series.E and self.rank == 6:
            return [
                [4/3,  5/3,  6/3,  4/3,  2/3,  3/3],
                [5/3, 10/3, 12/3,  8/3,  4/3,  6/3],
                [6/3, 12/3, 18/3, 12/3,  6/3,  9/3],
                [4/3,  8/3, 12/3, 10/3,  5/3,  6/3],
                [2/3,  4/3,  6/3,  5/3,  4/3,  3/3],
                [3/3,  6/3,  9/3,  6/3,  3/3,  6/3]
            ]

        if self.series == Series.E and self.rank == 7:
            return [
                [4/2,  6/2,  8/2,  6/2,  4/2,  2/2,  4/2],
                [6/2, 12/2, 16/2, 12/2,  8/2,  4/2,  8/2],
                [8/2, 16/2, 24/2, 18/2, 12/2,  6/2, 12/2],
                [6/2, 12/2, 18/2, 15/2, 10/2,  5/2,  9/2],
                [4/2,  8/2, 12/2, 10/2,  8/2,  4/2,  6/2],
                [2/2,  4/2,  6/2,  5/2,  4/2,  3/2,  3/2],
                [4/2,  8/2, 12/2,  9/2,  6/2,  3/2,  7/2]
            ]

        if self.series == Series.E and self.rank == 8:
            return [
                [4,   7, 10,  8,  6,  4,  2,  5],
                [7,  14, 20, 16, 12,  8,  4, 10],
                [10, 20, 30, 24, 18, 12,  6, 15],
                [8,  16, 24, 20, 15, 10,  5, 12],
                [6,  12, 18, 15, 12,  8,  4,  9],
                [4,   8, 12, 10,  8,  6,  3,  6],
                [2,   4,  6,  5,  4,  3,  2,  3],
                [5,  10, 15, 12,  9,  6,  3,  8]
            ]

        if self.series == Series.F:
            return [
                [2, 3,   2,   1],
                [3, 6,   4,   2],
                [2, 4,   3, 3/2],
                [1, 2, 3/2,   1]
            ]

        if self.series == Series.G:
            return [
                [2,   1],
                [1, 2/3]
            ]

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

        if self.series == Series.B:
            return (
                [
                    (2*self.rank - i) * (i + 1)
                    for i in range(self.rank - 1)
                ]
                + [int(self.rank * (self.rank + 1)/2)]
            )

        if self.series == Series.C:
            return [
                (2*self.rank - i - 1) * (i + 1)
                for i in range(self.rank)
            ]

        if self.series == Series.D:
            return (
                [
                    (2*self.rank - i - 2) * (i + 1)
                    for i in range(self.rank - 2)
                ]
                + [int(self.rank * (self.rank - 1)/2)] * 2
            )

        if self.series == Series.E:
            return {
                6: [16, 30, 42, 30, 16, 22],
                7: [34, 66, 96, 75, 52, 27, 49],
                8: [92, 182, 270, 220, 168, 114, 58, 136]
            }[self.rank]

        if self.series == Series.F:
            return [22, 42, 30, 16]

        if self.series == Series.G:
            return [10, 6]


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
        return Weight(itertools.chain.from_iterable(weights))

    @property
    @functools.lru_cache(maxsize=None)
    def level_vector(self):
        return Weight(
            itertools.chain.from_iterable(
                simple_algebra.level_vector
                for simple_algebra in self
            )
        )
