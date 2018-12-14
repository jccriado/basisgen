from invariants.algebras import Series, SimpleAlgebra
from invariants.representations import Irrep
from invariants.weights import Weight


algebra = SimpleAlgebra(Series.A, 2)

name_irrep_table = [
    ('1', Irrep(algebra, Weight([0, 0]))),
    ('3', Irrep(algebra, Weight([1, 0]))),
    ('6', Irrep(algebra, Weight([2, 0]))),
    ('8', Irrep(algebra, Weight([1, 1]))),
    ('10', Irrep(algebra, Weight([3, 0]))),
    ('15', Irrep(algebra, Weight([2, 1]))),
    ("15'", Irrep(algebra, Weight([4, 0]))),
    ('21', Irrep(algebra, Weight([0, 5]))),
    ('24', Irrep(algebra, Weight([1, 3]))),
    ('27', Irrep(algebra, Weight([2, 2])))
]


def SU3_irrep(n, m=None):
    if m is None:
        if n[-1] == '*':
            fst, snd = SU3_irrep(n[:-1]).highest_weight
            return SU3_irrep(snd, fst)
        else:
            return dict(name_irrep_table)[n]
    else:
        return Irrep(algebra, Weight([n, m]))


def SU3_name(irrep):
    try:
        return {irrep: name for name, irrep in name_irrep_table}[irrep]
    except KeyError:
        return {
            SU3_irrep(irrep.highest_weight[1], irrep.highest_weight[0]):
            name + "*"
            for name, irrep in name_irrep_table
        }[irrep]


def SU3_show(irreps):
    return " + ".join(
        f"({count}) {name}" if count > 1 else f"{name}"
        for name, count in (
                (SU3_name(irrep), count) for irrep, count in irreps.items()
        )
    )


singlet = SU3_irrep(0, 0)
triplet = SU3_irrep(1, 0)
anti_triplet = SU3_irrep(0, 1)
octet = SU3_irrep(1, 1)
sextet = SU3_irrep(2, 0)
anti_sextet = SU3_irrep(0, 2)
