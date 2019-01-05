import argparse
import cProfile

from basisgen.smeft import smeft, sm_field_classes


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
        '--include_lower_dimension',
        action='store_const',
        default=False,
        const=True
    )

    argument_parser.add_argument(
        '--detailed',
        action='store_const',
        default=False,
        const=True
    )

    argument_parser.add_argument(
        '--no_eom',
        action='store_const',
        default=False,
        const=True
    )

    argument_parser.add_argument(
        '--covariants',
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

    if arguments.covariants:
        operators_generator = smeft(arguments.number_of_flavors).covariants
    else:
        operators_generator = smeft(arguments.number_of_flavors).invariants

    operators = operators_generator(
        arguments.dimension,
        verbose=True,
        ignore_lower_dimension=not arguments.include_lower_dimension,
        use_eom=not arguments.no_eom
    )

    if arguments.profile:
        profiler.disable()

    if not arguments.covariants:
        print("Number of operators: {}".format(operators.count()))

    if arguments.detailed or arguments.covariants:
        print(operators)
    else:
        print(
            operators.show_by_classes(
                classes=sm_field_classes(arguments.number_of_flavors)
            )
        )

    if arguments.profile:
        profiler.print_stats(sort='time')
