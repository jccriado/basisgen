from invariants.fields import Operator, EFT
from invariants.smeft import sm_gauge_algebra, smeft, phi, phic, u, uc, gL, gR
from invariants.weights import Weight

import unittest
from collections import Counter


class TestSMEFT(unittest.TestCase):
    def test_dimension_5_basis(self):
        invariants_by_dimension = {2: 1, 4: 14, 5: 16}

        for dimension, number_of_invariants in invariants_by_dimension.items():
            self.assertEqual(
                smeft(1).invariants(dimension).count(),
                number_of_invariants
            )

    def test_higgs_dimension_8_basis(self):
        higgs_only = EFT(sm_gauge_algebra, [phi, phic])

        invariants = EFT.Invariants({
            Operator(Counter({phic: 1, phi: 1})): {0: 1},

            Operator(Counter({phic: 2, phi: 2})): {0: 1, 2: 2, 4: 3},

            Operator(Counter({phic: 3, phi: 3})): {0: 1, 2: 2},

            Operator(Counter({phic: 4, phi: 4})): {0: 1}
        })

        self.assertEqual(
            higgs_only.invariants(8, verbose=False),
            invariants
        )

    def test_3_flavors(self):
        self.assertEqual(smeft(3).invariants(4).count(), 62)

    def test_operator_u_uc_gL_gR_D(self):
        operator = u(1) * uc(1) * gL * gR

        self.assertEqual(
            operator.invariants(8, ignore_lower_dimensions=True),
            {1: 3}
        )

    def test_covariants_higgs(self):
        covariants = EFT(sm_gauge_algebra, [phi, phic]).covariants(3)

        known_covariants = {
            (Weight([0, 0, 0, 0, 1]), (-0.5,)):
            Counter({(phic._to_operator(), 0): 1, (phi * phic**2, 0): 1}),

            (Weight([1, 1, 0, 0, 1]), (-0.5,)):
            Counter({(phic._to_operator(), 1): 1})
        }

        for key in known_covariants:
            self.assertEqual(
                covariants[key],
                known_covariants[key]
            )


if __name__ == '__main__':
    unittest.main()
