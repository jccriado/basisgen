from invariants.algebras import Series, SimpleAlgebra, SemisimpleAlgebra
from invariants.weights import Weight

import unittest


class TestSimpleAlgebra(unittest.TestCase):
    def setUp(self):
        self.a1 = SimpleAlgebra(Series.A, 1)
        self.a2 = SimpleAlgebra(Series.A, 2)
        self.b4 = SimpleAlgebra(Series.B, 4)
        self.c6 = SimpleAlgebra(Series.C, 6)
        self.d5 = SimpleAlgebra(Series.D, 5)
        self.e6 = SimpleAlgebra(Series.E, 6)
        self.e7 = SimpleAlgebra(Series.E, 7)
        self.e8 = SimpleAlgebra(Series.E, 8)
        self.f4 = SimpleAlgebra(Series.F, 4)
        self.g2 = SimpleAlgebra(Series.G, 2)

    def test_cartan_matrices(self):
        known_cartan_matrices = {
            self.a1: [[2]],

            self.a2: [
                [2, -1],
                [-1, 2]
            ],

            self.b4: [
                [2, -1,  0,  0],
                [-1, 2, -1,  0],
                [0, -1,  2, -2],
                [0,  0, -1,  2]
            ],

            self.c6: [
                [2, -1,  0,  0,  0,  0],
                [-1, 2, -1,  0,  0,  0],
                [0, -1,  2, -1,  0,  0],
                [0,  0, -1,  2, -1,  0],
                [0,  0,  0, -1,  2, -1],
                [0,  0,  0,  0, -2,  2]
            ],

            self.d5: [
                [2, -1,  0,  0,  0],
                [-1, 2, -1,  0,  0],
                [0, -1,  2, -1, -1],
                [0,  0, -1,  2,  0],
                [0,  0, -1,  0,  2]
            ],

            self.e6: [
                [2, -1,  0,  0,  0,  0],
                [-1, 2, -1,  0,  0,  0],
                [0, -1,  2, -1,  0, -1],
                [0,  0, -1,  2, -1,  0],
                [0,  0,  0, -1,  2,  0],
                [0,  0, -1,  0,  0,  2]
            ],

            self.e7: [
                [2, -1,  0,  0,  0,  0,  0],
                [-1, 2, -1,  0,  0,  0,  0],
                [0, -1,  2, -1,  0,  0, -1],
                [0,  0, -1,  2, -1,  0,  0],
                [0,  0,  0, -1,  2, -1,  0],
                [0,  0,  0,  0, -1,  2,  0],
                [0,  0, -1,  0,  0,  0,  2]
            ],

            self.e8: [
                [2, -1,  0,  0,  0,  0,  0,  0],
                [-1, 2, -1,  0,  0,  0,  0,  0],
                [0, -1,  2, -1,  0,  0,  0, -1],
                [0,  0, -1,  2, -1,  0,  0,  0],
                [0,  0,  0, -1,  2, -1,  0,  0],
                [0,  0,  0,  0, -1,  2, -1,  0],
                [0,  0,  0,  0,  0, -1,  2,  0],
                [0,  0, -1,  0,  0,  0,  0,  2]
            ],

            self.f4: [
                [2, -1,  0,  0],
                [-1, 2, -2,  0],
                [0, -1,  2, -1],
                [0,  0, -1,  2]
            ],

            self.g2: [
                [2, -3],
                [-1, 2]
            ]
        }

        for algebra, cartan_matrix in known_cartan_matrices.items():
            self.assertEqual(algebra.cartan_matrix, cartan_matrix)

    def test_highest_root(self):
        known_highest_roots = {
            self.a1: Weight([2]),
            self.a2: Weight([1, 1]),
            self.b4: Weight([0, 1, 0, 0]),
            self.c6: Weight([2, 0, 0, 0, 0, 0]),
            self.d5: Weight([0, 1, 0, 0, 0]),
            self.e6: Weight([0, 0, 0, 0, 0, 1]),
            self.e7: Weight([1, 0, 0, 0, 0, 0, 0]),
            self.e8: Weight([0, 0, 0, 0, 0, 0, 1, 0]),
            self.f4: Weight([1, 0, 0, 0]),
            self.g2: Weight([1, 0])
        }

        for algebra, root in known_highest_roots.items():
            self.assertEqual(algebra.highest_root, root)

    def test_level_of_simple_roots(self):
        known_levels = {
            self.a1: 0,
            self.a2: 1,
            self.b4: 6,
            self.c6: 10,
            self.d5: 6,
            self.e6: 10,
            self.e7: 16,
            self.e8: 28,
            self.f4: 10,
            self.g2: 4
        }

        for algebra, level in known_levels.items():
            self.assertEqual(algebra.level_of_simple_roots, level)

    def test_positive_roots(self):
        known_positive_roots = {
            self.a2: {
                Weight([1, 1]),
                Weight([2, -1]), Weight([-1, 2])
            },

            SimpleAlgebra(Series.C, 2): {
                Weight([2, 0]),
                Weight([0, 1]),
                Weight([2, -1]), Weight([-2, 2])
            },

            SimpleAlgebra(Series.G, 2): {
                Weight([1, 0]),
                Weight([-1, 3]),
                Weight([0, 1]),
                Weight([1, -1]),
                Weight([2, -3]), Weight([-1, 2])
            },

            SimpleAlgebra(Series.A, 3): {
                Weight([1, 0, 1]),
                Weight([1, 1, -1]), Weight([-1, 1, 1]),
                Weight([2, -1, 0]), Weight([-1, 2, -1]), Weight([0, -1, 2])
            },

            SimpleAlgebra(Series.B, 3): {
                Weight([0, 1, 0]),
                Weight([1, -1, 2]),
                Weight([1, 0, 0]), Weight([-1, 0, 2]),
                Weight([1, 1, -2]), Weight([-1, 1, 0]),
                Weight([2, -1, 0]), Weight([-1, 2, -2]), Weight([0, -1, 2])
            },

            SimpleAlgebra(Series.C, 3): {
                Weight([2, 0, 0]),
                Weight([0, 1, 0]),
                Weight([1, -1, 1]), Weight([-2, 2, 0]),
                Weight([1, 1, -1]), Weight([-1, 0, 1]),
                Weight([2, -1, 0]), Weight([-1, 2, -1]), Weight([0, -2, 2])
            },

            SimpleAlgebra(Series.A, 4): {
                Weight([1, 0, 0, 1]),

                Weight([1, 0, 1, -1]), Weight([-1, 1, 0, 1]),

                Weight([1, 1, -1, 0]), Weight([-1, 1, 1, -1]),
                Weight([0, -1, 1, 1]),

                Weight([2, -1, 0, 0]), Weight([-1, 2, -1, 0]),
                Weight([0, -1, 2, -1]), Weight([0, 0, -1, 2])
            },

            SimpleAlgebra(Series.D, 5): {
                Weight([0, 1, 0, 0, 0]),

                Weight([1, -1, 1, 0, 0]),

                Weight([-1, 0, 1, 0, 0]), Weight([1, 0, -1, 1, 1]),

                Weight([-1, 1, -1, 1, 1]), Weight([1, 0, 0, -1, 1]),
                Weight([1, 0, 0, 1, -1]),

                Weight([0, -1, 0, 1, 1]), Weight([-1, 1, 0, -1, 1]),
                Weight([-1, 1, 0, 1, -1]), Weight([1, 0, 1, -1, -1]),

                Weight([0, -1, 1, -1, 1]), Weight([0, -1, 1, 1, -1]),
                Weight([-1, 1, 1, -1, -1]), Weight([1, 1, -1, 0, 0]),

                Weight([0, 0, -1, 0, 2]), Weight([0, 0, -1, 2, 0]),
                Weight([0, -1, 2, -1, -1]), Weight([-1, 2, -1, 0, 0]),
                Weight([2, -1, 0, 0, 0])
            }
        }

        for algebra, roots in known_positive_roots.items():
            self.assertEqual(set(algebra.positive_roots), roots)

    def test_metric(self):
        known_metrics = {
            self.a1: [[1/2]],

            self.a2: [
                [2/3, 1/3],
                [1/3, 2/3]
            ],

            self.b4: [
                [1,   1,   1, 1/2],
                [1,   2,   2,   1],
                [1,   2,   3, 3/2],
                [1/2, 1, 3/2,   1]
            ],

            self.c6: [
                [1/2, 1/2, 1/2, 1/2, 1/2, 1/2],
                [1/2,   1,   1,   1,   1,   1],
                [1/2,   1, 3/2, 3/2, 3/2, 3/2],
                [1/2,   1, 3/2,   2,   2,   2],
                [1/2,   1, 3/2,   2, 5/2, 5/2],
                [1/2,   1, 3/2,   2, 5/2,   3]
            ],

            self.d5: [
                [1,   1,   1, 1/2, 1/2],
                [1,   2,   2,   1,   1],
                [1,   2,   3, 3/2, 3/2],
                [1/2, 1, 3/2, 5/4, 3/4],
                [1/2, 1, 3/2, 3/4, 5/4]
            ],

            self.e6: [
                [4/3,  5/3,  6/3,  4/3,  2/3,  3/3],
                [5/3, 10/3, 12/3,  8/3,  4/3,  6/3],
                [6/3, 12/3, 18/3, 12/3,  6/3,  9/3],
                [4/3,  8/3, 12/3, 10/3,  5/3,  6/3],
                [2/3,  4/3,  6/3,  5/3,  4/3,  3/3],
                [3/3,  6/3,  9/3,  6/3,  3/3,  6/3]
            ],

            self.e7: [
                [2, 3,  4,    3, 2,   1,   2],
                [3, 6,  8,    6, 4,   2,   4],
                [4, 8, 12,    9, 6,   3,   6],
                [3, 6,  9, 15/2, 5, 5/2, 9/2],
                [2, 4,  6,    5, 4,   2,   3],
                [1, 2,  3,  5/2, 2, 3/2, 3/2],
                [2, 4,  6,  9/2, 3, 3/2, 7/2]
            ],

            self.e8: [
                [4,   7, 10,  8,  6,  4, 2,  5],
                [7,  14, 20, 16, 12,  8, 4, 10],
                [10, 20, 30, 24, 18, 12, 6, 15],
                [8,  16, 24, 20, 15, 10, 5, 12],
                [6,  12, 18, 15, 12,  8, 4,  9],
                [4,   8, 12, 10,  8,  6, 3,  6],
                [2,   4,  6,  5,  4,  3, 2,  3],
                [5,  10, 15, 12,  9,  6, 3,  8]
            ],

            self.f4: [
                [2, 3,   2,   1],
                [3, 6,   4,   2],
                [2, 4,   3, 3/2],
                [1, 2, 3/2,   1]
            ],

            self.g2: [
                [2,   1],
                [1, 2/3]
            ]
        }

        for algebra, metric in known_metrics.items():
            self.assertEqual(algebra.metric, metric)


class TestSemisimpleAlgebra(unittest.TestCase):
    def setUp(self):
        self.a1 = SimpleAlgebra(Series.A, 1)
        self.a2 = SimpleAlgebra(Series.A, 2)
        self.a9 = SimpleAlgebra(Series.A, 9)
        self.c7 = SimpleAlgebra(Series.C, 7)
        self.e7 = SimpleAlgebra(Series.E, 7)
        self.g2 = SimpleAlgebra(Series.G, 2)

        self.algebras = [self.a1, self.a2, self.a9, self.c7, self.e7, self.g2]
        self.algebra = SemisimpleAlgebra(self.algebras)

    def test_sum(self):
        self.assertEqual(
            sum(self.algebras, SemisimpleAlgebra([])),
            self.algebra
        )

    def test_split_weight(self):
        self.assertEqual(
            self.algebra.split_weight(Weight([
                1,
                2, 5,
                -8, 3, 4, 5, 4, 2, 11, 2, 4,
                0, -2, -3, -4, 5, 7, 6,
                4, 3, 4, 3, 4, 3, 4,
                0, 0
            ])),
            [
                Weight([1]),
                Weight([2, 5]),
                Weight([-8, 3, 4, 5, 4, 2, 11, 2, 4]),
                Weight([0, -2, -3, -4, 5, 7, 6]),
                Weight([4, 3, 4, 3, 4, 3, 4]),
                Weight([0, 0])
            ]
        )


if __name__ == '__main__':
    unittest.main()
