from invariants.algebras import Series, SimpleAlgebra, SemisimpleAlgebra
from invariants.parsing import parse_weight, parse_algebra
from invariants.weights import Weight

import unittest


class TestParsing(unittest.TestCase):
    def test_weights(self):
        self.assertEqual(
            parse_weight('3'),
            Weight([3])
        )

        self.assertEqual(
            parse_weight('1 2'),
            Weight([1, 2])
        )

        self.assertEqual(
            parse_weight('8 -9 1'),
            Weight([8, -9, 1])
        )

        self.assertEqual(
            parse_weight('  1  8 -5   -11  4 2 '),
            Weight([1, 8, -5, -11, 4, 2])
        )

    def test_algebras(self):
        self.assertEqual(
            parse_algebra('A1'),
            SimpleAlgebra(Series.A, 1)
        )

        self.assertEqual(
            parse_algebra('E7'),
            SimpleAlgebra(Series.E, 7)
        )

        self.assertEqual(
            parse_algebra('D123'),
            SimpleAlgebra(Series.D, 123)
        )

        self.assertEqual(
            parse_algebra('Sp8'),
            SimpleAlgebra(Series.C, 4)
        )

        self.assertEqual(
            parse_algebra('SO4 x SU3 x SU2'),
            SemisimpleAlgebra([
                SimpleAlgebra(Series.A, 1),
                SimpleAlgebra(Series.A, 1),
                SimpleAlgebra(Series.A, 2),
                SimpleAlgebra(Series.A, 1)
            ])
        )

        self.assertEqual(
            parse_algebra(' A3+ F4+B23 +  A3  '),
            SemisimpleAlgebra([
                SimpleAlgebra(Series.A, 3),
                SimpleAlgebra(Series.F, 4),
                SimpleAlgebra(Series.B, 23),
                SimpleAlgebra(Series.A, 3)
            ])
        )


if __name__ == '__main__':
    unittest.main()
