from invariants.algebras import Series, SimpleAlgebra
from invariants.representations import Irrep
from invariants.weights import Weight


algebra = SimpleAlgebra(Series.A, 1)


def irrep(n):
    return Irrep(algebra, Weight([n - 1]))


singlet = irrep(1)
doublet = irrep(2)
triplet = irrep(3)
