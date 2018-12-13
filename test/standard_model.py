import argparse
import cProfile

from invariants.fields import Field, EFT
from invariants.statistics import Statistics

import invariants.SU2 as SU2
import invariants.SU3 as SU3
from invariants.Lorentz import algebra as Lorentz_algebra
from invariants.Lorentz import scalar, L_spinor, R_spinor, L_tensor, R_tensor


SM_algebra = Lorentz_algebra + SU3.algebra + SU2.algebra

q = Field(
    'q', L_spinor, SU3.triplet + SU2.doublet, [1/6], Statistics.FERMION, 1.5
)
qc = Field(
    'qc', R_spinor, SU3.anti_triplet + SU2.doublet, [-1/6], Statistics.FERMION,
    1.5
)

u = Field(
    'u', R_spinor, SU3.triplet + SU2.singlet, [2/3], Statistics.FERMION, 1.5
)
uc = Field(
    'uc', L_spinor, SU3.anti_triplet + SU2.singlet, [-2/3], Statistics.FERMION,
    1.5
)

d = Field(
    'd', R_spinor, SU3.triplet + SU2.singlet, [-1/3], Statistics.FERMION, 1.5
)
dc = Field(
    'dc', L_spinor, SU3.anti_triplet + SU2.singlet, [1/3], Statistics.FERMION,
    1.5
)

l = Field(
    'l', L_spinor, SU3.singlet + SU2.doublet, [-1/2], Statistics.FERMION, 1.5
)
lc = Field(
    'lc', R_spinor, SU3.singlet + SU2.doublet, [1/2], Statistics.FERMION, 1.5
)

e = Field(
    'e', R_spinor, SU3.singlet + SU2.singlet, [-1], Statistics.FERMION, 1.5
)
ec = Field(
    'ec', L_spinor, SU3.singlet + SU2.singlet, [1], Statistics.FERMION, 1.5
)

phi = Field(
    'phi', scalar, SU3.singlet + SU2.doublet, [1/2], Statistics.BOSON, 1
)
phic = Field(
    'phic', scalar, SU3.singlet + SU2.doublet, [-1/2], Statistics.BOSON, 1
)

bL = Field(
    'bL', L_tensor, SU3.singlet + SU2.singlet, [0], Statistics.BOSON, 2
)
bR = Field(
    'bR', R_tensor, SU3.singlet + SU2.singlet, [0], Statistics.BOSON, 2
)

wL = Field(
    'wL', L_tensor, SU3.singlet + SU2.triplet, [0], Statistics.BOSON, 2
)
wR = Field(
    'wR', R_tensor, SU3.singlet + SU2.triplet, [0], Statistics.BOSON, 2
)

gL = Field(
    'gL', L_tensor, SU3.octet + SU2.singlet, [0], Statistics.BOSON, 2
)
gR = Field(
    'gR', R_tensor, SU3.octet + SU2.singlet, [0], Statistics.BOSON, 2
)

SM_fields = [
    q, qc, u, uc, d, dc, l, lc, e, ec,
    bL, bR, wL, wR, gL, gR, phi, phic
]

SMEFT = EFT(SM_algebra, SM_fields)

if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(
        description="Compute bases for the SMEFT"
    )

    argument_parser.add_argument(
        '--dimension',
        type=int,
        metavar='D',
        default=6,
        help='maximum dimension for the operators'
    )

    argument_parser.add_argument(
        '--profile',
        action='store_const',
        const=True,
        default=False
    )

    arguments = argument_parser.parse_args()

    if arguments.profile:
        profiler = cProfile.Profile()
        profiler.enable()
    
    invariants = SMEFT.invariants(arguments.dimension)

    if arguments.profile:
        profiler.disable()
    
    print("Number of invariants: {}".format(
        sum(
            count
            for counter in invariants.values()
            for count in counter.values()
        )
    ))

    print(EFT.show_invariants(invariants))

    if arguments.profile:
        profiler.print_stats(sort='time')
