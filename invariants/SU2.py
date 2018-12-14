from invariants.algebras import Series, SimpleAlgebra
from invariants.representations import Irrep
from invariants.weights import Weight


algebra = SimpleAlgebra(Series.A, 1)


def SU2_irrep(n):
    return Irrep(algebra, Weight([int(n) - 1]))


singlet = SU2_irrep(1)
doublet = SU2_irrep(2)
triplet = SU2_irrep(3)
