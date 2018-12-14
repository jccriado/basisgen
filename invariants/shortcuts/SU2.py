from invariants.algebras import Series, SimpleAlgebra
from invariants.representations import Irrep
from invariants.weights import Weight


SU2_algebra = SimpleAlgebra(Series.A, 1)


def SU2_irrep(n):
    return Irrep(SU2_algebra, Weight([int(n) - 1]))


SU2_singlet = SU2_irrep(1)
SU2_doublet = SU2_irrep(2)
SU2_triplet = SU2_irrep(3)
