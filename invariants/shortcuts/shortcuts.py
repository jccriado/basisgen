from invariants.representations import Irrep
from invariants.shortcuts.parsing import parse_algebra, parse_weight
from invariants.statistics import Statistics


algebra = parse_algebra


def irrep(algebra, highest_weight):
    return Irrep(parse_algebra(algebra), parse_weight(highest_weight))


boson = Statistics.BOSON
fermion = Statistics.FERMION
