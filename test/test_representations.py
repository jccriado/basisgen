from invariants.algebras import Series, SimpleAlgebra, SemisimpleAlgebra
from invariants.representations import Representation, Irrep
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
            self.assertEqual(irrep.representation.weights, weights)

    def test_su3_tensor_products(self):
        su3 = SimpleAlgebra(Series.A, 2)

        def irrep(name):
            if name[-1] == '*':
                conjugate = irrep(name[:-1])
                conjugate.highest_weight = Weight([
                    conjugate.highest_weight[1],
                    conjugate.highest_weight[0],
                ])
                return conjugate

            return {
                '1': Irrep(su3, Weight([0, 0])),
                '3': Irrep(su3, Weight([1, 0])),
                '6': Irrep(su3, Weight([2, 0])),
                '8': Irrep(su3, Weight([1, 1])),
                '10': Irrep(su3, Weight([3, 0])),
                '15': Irrep(su3, Weight([2, 1])),
                "15'": Irrep(su3, Weight([4, 0])),
                '21': Irrep(su3, Weight([0, 5])),
                '24': Irrep(su3, Weight([1, 3])),
                '27': Irrep(su3, Weight([2, 2]))
            }[name]

        known_decompositions = [
            (
                irrep('3*') * irrep('3*'),
                [irrep('3'), irrep('6*')]
            ),
            (
                irrep('3') * irrep('3*'),
                [irrep('1'), irrep('8')]
            ),
            (
                irrep('6') * irrep('3'),
                [irrep('8'), irrep('10')]
            ),
            (
                irrep('6') * irrep('3*'),
                [irrep('3'), irrep('15')]
            ),
            (
                irrep('6') * irrep('6'),
                [irrep('6*'), irrep('15'), irrep("15'")]
            ),
            (
                irrep('6') * irrep('6*'),
                [irrep('1'), irrep('8'), irrep('27')]
            ),
            (
                irrep('8') * irrep('3'),
                [irrep('3'), irrep('6*'), irrep('15')],
            ),
            (
                irrep('8') * irrep('6*'),
                [irrep('3'), irrep('6*'), irrep('15'), irrep('24')]
            ),
            (
                irrep('8') * irrep('8'),
                [irrep('1'), irrep('8'), irrep('8'),
                 irrep('10'), irrep('10*'), irrep('27')]
            )
        ]

        for product, decomposition in known_decompositions:
            self.assertEqual(product, collections.Counter(decomposition))


if __name__ == '__main__':
    unittest.main()
