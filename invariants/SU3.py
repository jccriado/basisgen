from invariants.algebras import Series, SimpleAlgebra
from invariants.representations import Irrep
from invariants.weights import Weight


algebra = SimpleAlgebra(Series.A, 2)


def irrep(n, m):
    return Irrep(algebra, Weight([n, m]))


singlet = irrep(0, 0)
triplet = irrep(1, 0)
anti_triplet = irrep(0, 1)
octet = irrep(1, 1)
sextet = irrep(2, 0)
anti_sextet = irrep(0, 2)
