import argparse
import cProfile

from invariants.fields import EFT
from invariants.smeft import (
    smeft, phi, phic, bL, bR, wL, wR, gL, gR, q, qc, L, Lc, u, uc, d, dc, e, ec
)


def parse_arguments():
    argument_parser = argparse.ArgumentParser(
        description="Compute bases for the SMEFT"
    )

    argument_parser.add_argument(
        '--dimension',
        type=int,
        metavar='d',
        default=6,
        help='maximum dimension for the operators'
    )

    argument_parser.add_argument(
        '--number_of_flavors',
        type=int,
        metavar='Nf',
        default=1,
        help='Number of fermion flavors'
    )

    argument_parser.add_argument(
        '--profile',
        action='store_const',
        const=True,
        default=False
    )

    argument_parser.add_argument(
        '--ignore_lower_dimension',
        action='store_const',
        default=False,
        const=True
    )

    return argument_parser.parse_args()


if __name__ == '__main__':
    arguments = parse_arguments()

    if arguments.profile:
        profiler = cProfile.Profile()
        profiler.enable()

    invariants = smeft(arguments.number_of_flavors).invariants(
        arguments.dimension,
        verbose=True,
        ignore_lower_dimension=arguments.ignore_lower_dimension
    )

    fermions = [
        fermion(1)
        for fermion in [q, qc, L, Lc, u, uc, d, dc, e, ec]
    ]

    gauge_bosons = [bL, bR, wL, wR, gL, gR]

    if arguments.profile:
        profiler.disable()

    print("Number of invariants: {}".format(EFT.count_invariants(invariants)))

    classes = {fermion: 'psi' for fermion in fermions}
    classes.update({gauge_boson: 'X' for gauge_boson in gauge_bosons})
    classes.update({phi: 'phi', phic: 'phi'})
    print(EFT.show_invariants(invariants, by_lines=False, classes=classes))

    if arguments.profile:
        profiler.print_stats(sort='time')
