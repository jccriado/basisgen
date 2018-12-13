from invariants.algebras import Series, SimpleAlgebra, SemisimpleAlgebra
from invariants.fields import Field
from invariants.statistics import Statistics
from invariants.parsing import (
    parse_weight, parse_algebra, parse_lorentz_highest_weight
)
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
            parse_algebra('Lorentz x SU3 x SU2'),
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

    def test_lorentz_highest_weight(self):
        self.assertEqual(
            parse_lorentz_highest_weight('L_spinor'),
            Weight([1, 0])
        )

        self.assertEqual(
            parse_lorentz_highest_weight('3 2'),
            Weight([3, 2])
        )

    def test_fields(self):
        phi_dict = {
            'lorentz_irrep': 'scalar',
            'internal_irrep': '0 0 1',
            'charges': [1/2],
            'statistics': 'boson',
            'dimension': 1
        }
        phi = Field.from_dict('phi', parse_algebra('SU3xSU2'), phi_dict)

        self.assertEqual(phi.name, 'phi')

        self.assertEqual(
            phi.internal_irrep.algebra,
            SemisimpleAlgebra([
                SimpleAlgebra(Series.A, 2),
                SimpleAlgebra(Series.A, 1)
            ])
        )

        self.assertEqual(phi.lorentz_irrep.highest_weight, Weight([0, 0]))

        self.assertEqual(phi.internal_irrep.highest_weight, Weight([0, 0, 1]))

        self.assertEqual(phi.charges, [1/2])

        self.assertEqual(phi.statistics, Statistics.BOSON)

        self.assertEqual(phi.dimension, 1)

        self.assertEqual(phi.number_of_derivatives, 0)


if __name__ == '__main__':
    unittest.main()
