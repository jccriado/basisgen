from invariants.representations import Irrep, IrrepCounter
from invariants.shortcuts import lorentz_algebra, vector
from invariants.statistics import Statistics
from invariants.partitions import partitions
from invariants.weights import Weight

from collections import Counter
import functools
import itertools
import math
from operator import mul


class Field(object):
    def __init__(
            self,
            name,
            lorentz_irrep,
            internal_irrep,
            charges=None,
            statistics=Statistics.BOSON,
            dimension=1,
            number_of_derivatives=0,
            number_of_flavors=1,
    ):
        if charges is None:
            charges = []

        self.name = name
        self.lorentz_irrep = lorentz_irrep
        self.internal_irrep = internal_irrep
        self.charges = charges
        self.statistics = statistics
        self.dimension = dimension
        self.number_of_derivatives = number_of_derivatives
        self.number_of_flavors = number_of_flavors

    def __hash__(self):
        return hash((self.name, self.number_of_derivatives))

    def __eq__(self, other):
        return (
            self.name == other.name
            and self.number_of_derivatives == other.number_of_derivatives
            and self.irrep == other.irrep
        )

    def __str__(self):
        if self.number_of_derivatives == 0:
            return self.name

        elif self.number_of_derivatives == 1:
            return "D{name}".format(name=self.name)

        else:
            return "D^{number_of_derivatives}({name})".format(
                name=self.name,
                number_of_derivatives=self.number_of_derivatives
            )

    __repr__ = __str__

    def __mul__(self, other):
        return self._to_operator() * other._to_operator()

    def __pow__(self, exponent):
        return Operator({self: exponent})

    def _to_operator(self):
        return Operator([self])

    def power_irreps(self, exponent):
        return (
            [
                self.irrep.power(inner_exponent, self.statistics)
                for inner_exponent in partition
            ]
            for partition in partitions(exponent, self.number_of_flavors)
        )

    @property
    def irrep(self):
        return self.lorentz_irrep + self.internal_irrep

    def differentiate(self, times, use_eom=True):
        if use_eom:
            highest_weight = (
                Weight([times, times]) + self.lorentz_irrep.highest_weight
            )
            lorentz_irreps = [Irrep(lorentz_algebra, highest_weight)]

        else:
            # TODO: check this
            lorentz_irreps = itertools.chain.from_iterable(
                self.lorentz_irrep
                * vector.power(n, statistics=Statistics.BOSON)
                for n in range(times % 2, times + 1, 2)
            )

        return {
            Field(
                name=self.name,
                lorentz_irrep=lorentz_irrep,
                internal_irrep=self.internal_irrep,
                charges=self.charges,
                statistics=self.statistics,
                dimension=self.dimension+times,
                number_of_derivatives=self.number_of_derivatives+times,
                number_of_flavors=self.number_of_flavors
            )
            for lorentz_irrep in lorentz_irreps
        }


class Operator(object):
    def __init__(self, content):
        self.content = Counter(content)

    def __hash__(self):
        return hash(frozenset(self.content.items()))

    def __eq__(self, other):
        return self.content == other.content

    def __str__(self):
        if not self.content:
            return "1"

        def power_str(field, exponent):
            if exponent == 1:
                return str(field)
            else:
                return "({})^{}".format(field, exponent)

        return " ".join(map(power_str, *zip(*self.content.items())))

    __repr__ = __str__

    def __mul__(self, other):
        return Operator(self.content + other._to_operator().content)

    def _to_operator(self):
        return self

    @property
    def dimension(self):
        return sum(
            field.dimension * exponent
            for field, exponent in self.content.items()
        )

    @property
    def charges(self):
        transposed_charges = zip(*(
            [charge * exponent for charge in field.charges]
            for field, exponent in self.content.items()
        ))

        return [sum(charges) for charges in transposed_charges]

    @property
    def is_neutral(self):
        return all(
            math.isclose(charge, 0, abs_tol=1e-10)
            for charge in self.charges
        )

    def differentiate_fields(self, times):
        def differentiate_field_by_partition(field, partition):
            chains = itertools.product(*(
                field.differentiate(number_of_derivatives)
                for number_of_derivatives in partition
            ))

            return (Operator(chain) for chain in chains)

        def chain_products(iterable_of_iterables):
            return itertools.chain.from_iterable(
                itertools.product(*inner_iterable)
                for inner_iterable in iterable_of_iterables
            )

        diff_partitions = chain_products(
            (
                partitions(number_of_derivatives, exponent)
                for exponent, number_of_derivatives in
                zip(self.content.values(), content_partition)
            )
            for content_partition in partitions(times, len(self.content))
        )

        operators_by_partition = chain_products(
            (
                differentiate_field_by_partition(field, partition)
                for partition, field in zip(content_partition, self.content)
            )
            for content_partition in set(diff_partitions)
        )

        return set(
            functools.reduce(mul, operators)
            for operators in operators_by_partition
        )

    @staticmethod
    def total_derivatives(initial_irrep, max_derivatives, initial_derivatives):
        def possible_irreps(lorentz_irrep):
            return initial_irrep * (
                lorentz_irrep + Irrep.singlet(initial_irrep.algebra[2:])
            )

        def vector_power(number_of_derivatives):
            return vector.power(number_of_derivatives, Statistics.BOSON)

        return {
            initial_derivatives + number_of_derivatives:
            IrrepCounter(
                irrep
                for lorentz_irrep in vector_power(number_of_derivatives)
                for irrep in possible_irreps(lorentz_irrep)
            )
            for number_of_derivatives in range(1, max_derivatives + 1)
        }

    @property
    def irreps(self):
        possibilities = itertools.product(*(
            field.power_irreps(exponent)
            for field, exponent in self.content.items()
        ))

        chains = map(itertools.chain.from_iterable, possibilities)

        return sum(
            (functools.reduce(mul, chain) for chain in chains),
            Counter()
        )

    def irreps_with_derivatives(self, max_dimension, filter_internal_singlets):
        max_derivatives = int(max_dimension - self.dimension)

        return {
            n_derivatives:
            IrrepCounter(
                irrep
                for operator in self.differentiate_fields(n_derivatives)
                for irrep in operator.irreps.elements()
                if not filter_internal_singlets or irrep[2:].is_singlet
            )
            for n_derivatives in range(max_derivatives + 1)
        }

    def irreps_without_total_derivatives(self, max_dimension):
        def total_derivatives(initial_irrep, derivative_count):
            return Operator.total_derivatives(
                initial_irrep,
                int(max_dimension - self.dimension - derivative_count),
                derivative_count
            )

        def remove_tower(irreps_by_derivatives, tower, initial_count):
            for n_derivatives, counter in tower.items():
                irreps_by_derivatives[n_derivatives] -= IrrepCounter({
                    irrep: count * initial_count
                    for irrep, count in counter.items()
                })

        out_irreps = self.irreps_with_derivatives(max_dimension, True)

        for derivative_count in range(int(max_dimension - self.dimension)):
            current_irreps = out_irreps[derivative_count].items()

            for initial_irrep, initial_count in current_irreps:
                tower = total_derivatives(initial_irrep, derivative_count)
                remove_tower(out_irreps, tower, initial_count)

        return out_irreps

    def invariants(self, max_dimension, ignore_lower_dimensions=False):
        def correct_dimension(number_of_derivatives):
            return (
                not ignore_lower_dimensions
                or number_of_derivatives == max_dimension - self.dimension
            )

        def sum_singlets(irrep_counter):
            return sum(
                count for irrep, count in irrep_counter.items()
                if irrep.is_singlet
            )

        irreps = self.irreps_without_total_derivatives(max_dimension)

        return {
            number_of_derivatives: sum_singlets(irrep_counter)
            for number_of_derivatives, irrep_counter in irreps.items()
            if correct_dimension(number_of_derivatives)
            if sum_singlets(irrep_counter)
        }

    def covariants(self, max_dimension, ignore_lower_dimensions=False):
        irreps = self.irreps_with_derivatives(max_dimension, False)
        max_derivatives = max_dimension - self.dimension

        return {
            n_derivatives: current_irreps
            for n_derivatives, current_irreps in irreps.items()
            if not ignore_lower_dimensions or n_derivatives == max_derivatives
            if irreps
        }


class EFT(object):
    class Invariants(object):
        def __init__(self, invariants):
            self.invariants = invariants

        def __eq__(self, other):
            return self.invariants == other.invariants

        @staticmethod
        def _show_derivatives(n_derivatives):
            if n_derivatives == 0:
                return ""
            elif n_derivatives == 1:
                return " D"
            else:
                return " D^" + str(n_derivatives)

        @staticmethod
        def _show_item(count, operator, n_derivatives):
            return "{operator}{derivatives}: {count}".format(
                count=count,
                operator=operator,
                derivatives=EFT.Invariants._show_derivatives(n_derivatives)
            )

        def __str__(self):
            return ("\n").join(
                EFT.Invariants._show_item(count, op, n_derivatives)
                for op, current_invariants in self.invariants.items()
                for n_derivatives, count in current_invariants.items()
            )

        def show_by_classes(self, classes):
            def power_str(field, exponent):
                if exponent == 1:
                    return str(field)
                else:
                    return "({})^{}".format(field, exponent)

            contents = {}
            for operator, ders in self.invariants.items():
                content_chain = (
                    Counter({classes[field]: exponent})
                    for field, exponent in operator.content.items()
                )

                content = sum(content_chain, Counter()).items()

                content_str = " ".join(
                    power_str(field, exponent) for field, exponent in content
                )

                contents.setdefault(content_str, Counter())
                contents[content_str] += ders

            return "\n".join(
                EFT.Invariants._show_item(count, op, n_derivatives)
                for op, current_invariants in contents.items()
                for n_derivatives, count in current_invariants.items()
            )

        def count(self):
            return sum(
                count
                for counter in self.invariants.values()
                for count in counter.values()
            )

    class Covariants(object):
        def __init__(self, covariants):
            self.covariants = covariants

        def __eq__(self, other):
            return self.covariants == other.covariants

        @staticmethod
        def _show_item(count, operator, n_derivatives):
            return EFT.Invariants._show_item(
                count, operator, n_derivatives
            )

        @staticmethod
        def _show_irrep_charges(irrep, charges):
            if len(charges) == 1:
                return "[irrep={}, charge={}]".format(irrep, charges[0])
            else:
                return "[irrep={}, charges={}]".format(irrep, charges)

        @staticmethod
        def _show_irrep_operators(irrep_charges, operator_counter):
            return "{irrep}: {operators}".format(
                irrep=EFT.Covariants._show_irrep_charges(*irrep_charges),
                operators=" + ".join(
                    EFT.Covariants._show_item(count, op, n_derivatives)
                    for (op, n_derivatives), count in operator_counter.items()
                )
            )

        def __str__(self):
            return "\n".join(
                EFT.Covariants._show_irrep_operators(key, operator_counter)
                for key, operator_counter in self.covariants.items()
            )

        def __getitem__(self, key):
            return self.covariants[key]

    def __init__(self, internal_algebra, fields, use_eom=False):
        self.algebra = lorentz_algebra + internal_algebra
        self.fields = fields
        self.use_eom = use_eom
        self.cached_results = {}

    @staticmethod
    def _combinations(fields, max_dimension):
        if not fields:
            return [Counter()]

        max_exponent = math.floor(max_dimension / fields[0].dimension)

        return [
            Counter({fields[0]: exponent}) + combination
            for exponent in range(max_exponent + 1)
            for combination in EFT._combinations(
                    fields[1:],
                    max_dimension - exponent * fields[0].dimension
            )
        ]

    def operators(self, max_dimension):
        return map(Operator, EFT._combinations(self.fields, max_dimension))

    def invariants(
            self,
            max_dimension,
            verbose=False,
            ignore_lower_dimension=False
    ):
        result = {}

        operators_printer = ProgressPrinter(
            "Computing field content combinations...", "done.",
            printing_function=print if verbose else None
        )
        invariants_printer = ProgressPrinter(
            "Computing invariants...", "done.", "({progress}/{total})",
            printing_function=print if verbose else None
        )

        operators_printer.start()
        operators = list(self.operators(max_dimension))
        total = len(operators)
        operators_printer.end()

        invariants_printer.start()
        for progress, operator in enumerate(operators):
            invariants_printer.update(progress=progress, total=total)

            if operator.content and operator.is_neutral:
                result[operator] = operator.invariants(
                    max_dimension,
                    ignore_lower_dimension
                )
        invariants_printer.end()

        return EFT.Invariants(result)

    def covariants(
            self,
            max_dimension,
            verbose=False,
            ignore_lower_dimension=False
    ):
        result = {}

        operators_printer = ProgressPrinter(
            "Computing field content combinations...", "done.",
            printing_function=print if verbose else None
        )
        covariants_printer = ProgressPrinter(
            "Computing covariant operators...", "done.",
            "({progress}/{total})",
            printing_function=print if verbose else None
        )

        operators_printer.start()
        operators = list(self.operators(max_dimension))
        total = len(operators)
        operators_printer.end()

        covariants_printer.start()
        for progress, operator in enumerate(operators):
            covariants_printer.update(progress=progress, total=total)

            if operator.content:
                covariants = operator.covariants(
                    max_dimension,
                    ignore_lower_dimension
                )
                for number_of_derivatives, irreps in covariants.items():
                    for irrep, count in irreps.items():
                        irrep_with_charges = (
                            irrep.highest_weight,
                            tuple(operator.charges)
                        )
                        result.setdefault(irrep_with_charges, Counter())
                        result[irrep_with_charges] += Counter({
                            (operator, number_of_derivatives): count
                        })
        covariants_printer.end()

        return EFT.Covariants(result)


class ProgressPrinter(object):
    def __init__(
            self,
            starting_message,
            ending_message,
            progress_template=None,
            printing_function=print
    ):
        self.starting_message = starting_message
        self.ending_message = ending_message
        self.progress_template = progress_template
        self.max_message_length = len(self.starting_message)

        if printing_function is None:
            def do_nothing(*args, **kwargs):
                return

            self.printing_function = do_nothing

        else:
            self.printing_function = printing_function

        self.start()

    def start(self):
        self.printing_function(self.starting_message, end='\r', flush=True)

    def end(self):
        self.clear()
        self.printing_function(
            self.starting_message + " " + self.ending_message,
            end='\n',
            flush=True
        )

    def update(self, **kwargs):
        message = (
            self.starting_message + " " +
            self.progress_template.format(**kwargs)
        )

        self.max_message_length = len(message)
        self.printing_function(message, end='\r', flush=True)

    def clear(self):
        self.printing_function(
            " " * self.max_message_length,
            end='\r',
            flush=True
        )
