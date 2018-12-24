from basisgen.representations import Irrep
from basisgen.parsing import parse_algebra, parse_weight
from basisgen.statistics import Statistics


algebra = parse_algebra


def irrep(algebra, highest_weight):
    return Irrep(parse_algebra(algebra), parse_weight(highest_weight))


boson = Statistics.BOSON
fermion = Statistics.FERMION
