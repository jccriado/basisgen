from invariants.algebras import Series, SimpleAlgebra
from invariants.representations import Irrep
from invariants.weights import Weight

import unittest
import collections


class TestIrrep(unittest.TestCase):
    def test_representation(self):
        known_weights = {
            Irrep(
                SimpleAlgebra(Series.D, 5),
                Weight([0, 0, 0, 0, 1])
            ):
            collections.Counter([
                Weight([0, 0, 0, 0, 1]),
                Weight([0, 0, 1, 0, -1]),
                Weight([0, 1, -1, 1, 0]),
                Weight([1, -1, 0, 1, 0]), Weight([0, 1, 0, -1, 0]),
                Weight([-1, 0, 0, 1, 0]), Weight([1, -1, 1, -1, 0]),
                Weight([-1, 0, 1, -1, 0]), Weight([1, 0, -1, 0, 1]),
                Weight([-1, 1, -1, 0, 1]), Weight([1, 0, 0, 0, -1]),
                Weight([0, -1, 0, 0, 1]), Weight([-1, 1, 0, 0, -1]),
                Weight([0, -1, 1, 0, -1]),
                Weight([0, 0, -1, 1, 0]),
                Weight([0, 0, 0, -1, 0])
            ]),

            Irrep(
                SimpleAlgebra(Series.E, 6),
                Weight([1, 0, 0, 0, 0, 0])
            ):
            collections.Counter([
                Weight([1, 0, 0, 0, 0, 0]),
                Weight([-1, 1, 0, 0, 0, 0]),
                Weight([0, -1, 1, 0, 0, 0]),
                Weight([0, 0, -1, 1, 0, 1]),
                Weight([0, 0, 0, -1, 1, 1]), Weight([0, 0, 0, 1, 0, -1]),
                Weight([0, 0, 0, 0, -1, 1]), Weight([0, 0, 1, -1, 1, -1]),
                Weight([0, 0, 1, 0, -1, -1]), Weight([0, 1, -1, 0, 1, 0]),
                Weight([0, 1, -1, 1, -1, 0]), Weight([1, -1, 0, 0, 1, 0]),

                Weight([0, 1, 0, -1, 0, 0]), Weight([1, -1, 0, 1, -1, 0]),
                Weight([-1, 0, 0, 0, 1, 0]),

                Weight([1, -1, 1, -1, 0, 0]), Weight([-1, 0, 0, 1, -1, 0]),
                Weight([1, 0, -1, 0, 0, 1]), Weight([-1, 0, 1, -1, 0, 0]),
                Weight([1, 0, 0, 0, 0, -1]), Weight([-1, 1, -1, 0, 0, 1]),
                Weight([-1, 1, 0, 0, 0, -1]), Weight([0, -1, 0, 0, 0, 1]),
                Weight([0, -1, 1, 0, 0, -1]),
                Weight([0, 0, -1, 1, 0, 0]),
                Weight([0, 0, 0, -1, 1, 0]),
                Weight([0, 0, 0, 0, -1, 0])
            ])
        }

        for irrep, weights in known_weights.items():
            self.assertEqual(irrep.weight_system.weights, weights)

    def test_su3_tensor_products(self):
        SU3_algebra = SimpleAlgebra(Series.A, 2)

        name_irrep_table = {
            '1': Irrep(SU3_algebra, Weight([0, 0])),
            '3': Irrep(SU3_algebra, Weight([1, 0])),
            '6': Irrep(SU3_algebra, Weight([2, 0])),
            '8': Irrep(SU3_algebra, Weight([1, 1])),
            '10': Irrep(SU3_algebra, Weight([3, 0])),
            '15': Irrep(SU3_algebra, Weight([2, 1])),
            "15'": Irrep(SU3_algebra, Weight([4, 0])),
            '21': Irrep(SU3_algebra, Weight([0, 5])),
            '24': Irrep(SU3_algebra, Weight([1, 3])),
            '27': Irrep(SU3_algebra, Weight([2, 2]))
        }

        def SU3_irrep(name):
            if name in name_irrep_table:
                return name_irrep_table[name]
            else:
                return name_irrep_table[name[:-1]].conjugate

        known_decompositions = [
            (
                SU3_irrep('3*') * SU3_irrep('3*'),
                [SU3_irrep('3'), SU3_irrep('6*')]
            ),
            (
                SU3_irrep('3') * SU3_irrep('3*'),
                [SU3_irrep('1'), SU3_irrep('8')]
            ),
            (
                SU3_irrep('6') * SU3_irrep('3'),
                [SU3_irrep('8'), SU3_irrep('10')]
            ),
            (
                SU3_irrep('6') * SU3_irrep('3*'),
                [SU3_irrep('3'), SU3_irrep('15')]
            ),
            (
                SU3_irrep('6') * SU3_irrep('6'),
                [SU3_irrep('6*'), SU3_irrep('15'), SU3_irrep("15'")]
            ),
            (
                SU3_irrep('6') * SU3_irrep('6*'),
                [SU3_irrep('1'), SU3_irrep('8'), SU3_irrep('27')]
            ),
            (
                SU3_irrep('8') * SU3_irrep('3'),
                [SU3_irrep('3'), SU3_irrep('6*'), SU3_irrep('15')],
            ),
            (
                SU3_irrep('8') * SU3_irrep('6*'),
                [SU3_irrep('3'), SU3_irrep('6*'),
                 SU3_irrep('15'), SU3_irrep('24')]
            ),
            (
                SU3_irrep('8') * SU3_irrep('8'),
                [SU3_irrep('1'), SU3_irrep('8'), SU3_irrep('8'),
                 SU3_irrep('10'), SU3_irrep('10*'), SU3_irrep('27')]
            )
        ]

        for product, decomposition in known_decompositions:
            self.assertEqual(product, collections.Counter(decomposition))

    def test_positive_roots(self):
        known_positive_roots = {
            SimpleAlgebra(Series.A, 2): {
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
            self.assertEqual(set(Irrep.positive_roots(algebra)), roots)


if __name__ == '__main__':
    unittest.main()
