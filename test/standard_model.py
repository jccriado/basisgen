import argparse
import cProfile

from invariants.fields import EFT
from invariants.smeft import smeft

if __name__ == '__main__':
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
        '--profile',
        action='store_const',
        const=True,
        default=False
    )

    argument_parser.add_argument(
        '--number_of_flavors',
        type=int,
        metavar='Nf',
        default=1,
        help='Number of different fermion flavors'
    )

    arguments = argument_parser.parse_args()

    if arguments.profile:
        profiler = cProfile.Profile()
        profiler.enable()

    invariants = smeft(
        arguments.number_of_flavors
    ).invariants(arguments.dimension, verbose=True)

    if arguments.profile:
        profiler.disable()

    print("Number of invariants: {}".format(EFT.count_invariants(invariants)))

    print(EFT.show_invariants(invariants))

    if arguments.profile:
        profiler.print_stats(sort='time')
