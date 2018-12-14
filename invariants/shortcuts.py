from invariants.representations import Irrep
from invariants.parsing import parse_algebra, parse_weight

algebra = parse_algebra


def irrep(algebra, highest_weight):
    return Irrep(parse_algebra(algebra), parse_weight(highest_weight))
