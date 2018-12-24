from basisgen.algebras import Series, SimpleAlgebra, SemisimpleAlgebra
from basisgen.lorentz import lorentz_algebra
from basisgen.weights import Weight

import re


def parse_weight(code):
    return Weight(map(int, code.split()))


def _parse_simple_group(code):
    match = re.match(r'(?P<series>SU|SO|Sp)(?P<N>[0-9]+)', code)
    series = match.group('series')
    N = int(match.group('N'))

    if series == 'SU':
        return SimpleAlgebra(Series.A, N - 1)
    elif series == 'SO' and N % 2 == 1:
        return SimpleAlgebra(Series.B, round((N - 1) / 2))
    elif series == 'Sp':
        return SimpleAlgebra(Series.C, round(N/2))
    elif series == 'SO' and N == 4:
        return lorentz_algebra
    elif series == 'SO' and N % 2 == 0:
        return SimpleAlgebra(Series.D, round(N/2))


def _parse_simple_algebra(code):
    match = re.match(r'(?P<series>[ABCDEFG])(?P<rank>[0-9]+)', code)

    series = {
        'A': Series.A,
        'B': Series.B,
        'C': Series.C,
        'D': Series.D,
        'E': Series.E,
        'F': Series.F,
        'G': Series.G
    }[match.group('series')]

    return SimpleAlgebra(series, int(match.group('rank')))


def parse_algebra(code):
    code = code.replace(' ', '')

    if '+' in code and 'x' in code:
        raise Exception(
            "Mixed algebra-style and group-style notation in '{}'".format(code)
        )

    elif '+' in code:
        return sum(
            map(_parse_simple_algebra, code.split('+')),
            SemisimpleAlgebra([])
        )

    elif 'x' in code:
        return sum(
            map(_parse_simple_group, code.split('x')),
            SemisimpleAlgebra([])
        )

    else:
        try:
            return _parse_simple_algebra(code)
        except Exception:
            return _parse_simple_group(code)
