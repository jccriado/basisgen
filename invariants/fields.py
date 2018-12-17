from invariants.algebras import Series, SimpleAlgebra, SemisimpleAlgebra
from invariants.representations import Irrep
from invariants.shortcuts import lorentz_algebra
from invariants.statistics import Statistics
from invariants.partitions import partitions
from invariants.weights import Weight

from collections import Counter
import copy
import functools
import itertools
import math


_lorentz_algebra = SemisimpleAlgebra([
    SimpleAlgebra(Series.A, 1),
    SimpleAlgebra(Series.A, 1)
])

_lorentz_vector_irrep = Irrep(_lorentz_algebra, Weight([1, 1]))


class Field(object):
    def __init__(
            self,
            name,
            lorentz_irrep,
            internal_irrep,
            charges,
            statistics,
            dimension,
            number_of_derivatives=0,
            number_of_flavors=1,
    ):
        self.name = name
        self.lorentz_irrep = lorentz_irrep
        self.internal_irrep = internal_irrep
        self.charges = charges
        self.statistics = statistics
        self.dimension = dimension
        self.number_of_derivatives = number_of_derivatives
        self.number_of_flavors = number_of_flavors

    def __hash__(self):
        return hash(self.name)

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
            return f"D{self.name}"
        else:
            return f"D^{self.number_of_derivatives}({self.name})"

    __repr__ = __str__

    def __mul__(self, other):
        return self._to_operator() * other._to_operator()

    def __pow__(self, exponent):
        return Operator(Counter({self: exponent}))

    def _to_operator(self):
        return Operator(Counter([self]))

    @property
    def irrep(self):
        return self.lorentz_irrep + self.internal_irrep

    def differentiate(self, times, use_eom=False):
        if use_eom:
            highest_weight = (
                Weight([times, times])
                + self.lorentz_irrep.highest_weight
            )
            lorentz_irreps = [Irrep(_lorentz_algebra, highest_weight)]

        else:
            # TODO: check this
            lorentz_irreps = itertools.chain.from_iterable(
                Irrep(_lorentz_algebra, Weight([n, n])) * self.lorentz_irrep
                for n in range(times % 2, times + 1, 2)
            )

        return set(
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
        )


class Operator(object):
    def __init__(self, content):
        self.content = content

    def __hash__(self):
        return hash(frozenset(self.content.items()))

    def __eq__(self, other):
        return self.content == other.content

    def __str__(self):
        if not self.content:
            return "1"

        def power_str(item):
            field, exponent = item
            if exponent == 1:
                return str(field)
            else:
                return f"({field})^{exponent}"

        return " ".join(map(power_str, self.content.items()))

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
        return [
            sum(charges) for charges in zip(*[
                [
                    charge * exponent
                    for charge in field.charges
                ]
                for field, exponent in self.content.items()
            ])
        ]

    @property
    def is_neutral(self):
        # TODO: fix rounding errors
        return all(
            math.isclose(charge, 0, abs_tol=1e-05)
            for charge in self.charges
        )

    @property
    def irreps(self):
        power_irreps = list(
            [
                [
                    field.irrep.power(inner_exponent, field.statistics)
                    for inner_exponent in partition
                ]
                for partition in partitions(exponent, field.number_of_flavors)
            ]
            for field, exponent in self.content.items()
        )

        if not self.content:
            return Counter()

        return sum(
            (
                functools.reduce(
                    Irrep.product,
                    itertools.chain.from_iterable(power_irreps_instance)
                )
                for power_irreps_instance in itertools.product(*power_irreps)
            ),
            Counter()
        )

    @property
    def derivatives(self):
        results = set()
        for field in self.content:
            for new_field in field.differentiate(1, use_eom=True):
                new_content = copy.copy(self.content)
                new_content -= Counter({field: 1})
                new_content[new_field] += 1

                results.add(Operator(new_content))

        return results

    def differentiate(self, times):
        if times == 0:
            return {self}
        else:
            return {
                operator
                for derivative in self.derivatives
                for operator in derivative.differentiate(times - 1)
            }

    def internal_singlets(self, max_dimension, filter_internal_singlets):
        max_derivatives = int(max_dimension - self.dimension)

        return {
            number_of_derivatives: Counter(
                irrep
                for operator in self.differentiate(number_of_derivatives)
                for irrep in operator.irreps.elements()
                if not filter_internal_singlets or irrep[2:].is_singlet
            )
            for number_of_derivatives in range(max_derivatives + 1)
        }

    @staticmethod
    def descendants(initial_irrep, max_derivatives, initial_derivatives):
        return {
            initial_derivatives + number_of_derivatives: Counter(
                irrep
                for lorentz_irrep in _lorentz_vector_irrep.power(
                            number_of_derivatives,
                            Statistics.BOSON
                )
                for irrep in (
                        (lorentz_irrep
                         + Irrep.singlet(initial_irrep.algebra[2:]))
                        * initial_irrep
                )
            )
            for number_of_derivatives in range(1, max_derivatives + 1)
        }

    def irreps_by_derivatives(self, max_dimension):
        out_irreps = self.internal_singlets(
            max_dimension,
            filter_internal_singlets=True
        )

        for derivative_count in range(int(max_dimension - self.dimension)):
            current_irreps = out_irreps[derivative_count].items()
            for initial_irrep, initial_count in current_irreps:
                descendants = Operator.descendants(
                    initial_irrep,
                    int(max_dimension - self.dimension - derivative_count),
                    int(derivative_count)
                )

                for number_of_derivatives, counter in descendants.items():
                    out_irreps[number_of_derivatives] -= Counter({
                        irrep: count * initial_count
                        for irrep, count in counter.items()
                    })

        return out_irreps

    def invariants(self, max_dimension, ignore_lower_dimensions=False):
        irreps_by_derivatives = self.irreps_by_derivatives(max_dimension)
        result = {}
        for number_of_derivatives, counter in irreps_by_derivatives.items():
            if (
                    not ignore_lower_dimensions
                    or number_of_derivatives == max_dimension - self.dimension
            ):
                total = sum(
                    count
                    for irrep, count in counter.items() if irrep.is_singlet
                )
                if total > 0:
                    result[number_of_derivatives] = total

        return result

    def covariants(self, max_dimension, ignore_lower_dimensions=False):
        irreps_by_derivatives = self.internal_singlets(
            max_dimension,
            False
        )

        if ignore_lower_dimensions:
            number_of_derivatives = max_dimension - self.dimension
            irreps = irreps_by_derivatives[number_of_derivatives]
            if irreps:
                return {number_of_derivatives: irreps}
            else:
                return {}
        else:
            return {
                number_of_derivatives: irreps
                for number_of_derivatives, irreps
                in irreps_by_derivatives.items()
                if irreps
            }


class EFT(object):
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

    def _operators(self, max_dimension):
        return map(Operator, EFT._combinations(self.fields, max_dimension))

    def operators(self, max_dimension):
        return self.cached_results.setdefault(
            ('operators', max_dimension),
            self._operators(max_dimension)
        )

    def _invariants(self, max_dimension, verbose, ignore_lower_dimension):
        result = {}

        if verbose:
            print(
                "Computing field content combinations... ", end="", flush=True
            )

        operators = list(self.operators(max_dimension))
        total = len(operators)
        spaces = " " * (round(math.log(total, 10)) + 1)

        if verbose:
            print("done.")
            print(f"Computing invariants... (0/{total})", end="\r", flush=True)

        for progress, operator in enumerate(operators):
            if verbose and progress % round(total / 50) == 0:
                print(
                    f"Computing invariants ({progress}/{total})" + spaces,
                    end="\r",
                    flush=True
                )

            if operator.content and operator.is_neutral:
                result[operator] = operator.invariants(
                    max_dimension,
                    ignore_lower_dimension
                )

        if verbose:
            print(f"Computing invariants... done." + spaces * 2, flush=True)

        return result

    def invariants(
            self,
            max_dimension,
            verbose=False,
            ignore_lower_dimension=False
    ):
        return self.cached_results.setdefault(
            ('invariants', max_dimension),
            self._invariants(max_dimension, verbose, ignore_lower_dimension)
        )

    @staticmethod
    def show_invariants(invariants, by_lines=False, classes=None):
        if classes is None:
            contents = invariants
        else:
            def power_str(item):
                field, exponent = item
                if exponent == 1:
                    return str(field)
                else:
                    return f"({field})^{exponent}"

            contents = {}
            for operator, ders in invariants.items():
                content = " ".join(map(
                    power_str,
                    sum(
                            (
                                Counter({classes[field]: exponent})
                                for field, exponent in operator.content.items()
                            ),
                            Counter()
                    ).items()
                ))
                contents.setdefault(content, Counter())
                contents[content] += ders

        return ("\n" if by_lines else " + ").join(
            "{count}{operator}{derivatives}".format(
                count=str(count) + " " if count > 1 or by_lines else "",
                operator=operator,
                derivatives="" if number_of_derivatives == 0 else (
                    " D" if number_of_derivatives == 1 else
                    f" D^{number_of_derivatives}"
                )
            )
            for operator, current_invariants in contents.items()
            for number_of_derivatives, count in current_invariants.items()
        )

    @staticmethod
    def count_invariants(invariants, verbose=False):
        return sum(
            count
            for counter in invariants.values()
            for count in counter.values()
        )

    def _covariants(self, max_dimension, verbose, ignore_lower_dimension):
        result = {}

        if verbose:
            print(
                "Computing field content combinations... ", end="", flush=True
            )

        operators = list(self.operators(max_dimension))
        total = len(operators)
        spaces = " " * (round(math.log(total, 10)) + 1)

        if verbose:
            print("done.")
            print(f"Computing covariants... (0/{total})", end="\r", flush=True)

        for progress, operator in enumerate(operators):
            if verbose and progress % round(total / 50) == 0:
                print(
                    f"Computing covariants ({progress}/{total})" + spaces,
                    end="\r",
                    flush=True
                )

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

        if verbose:
            print(f"Computing covariants... done." + spaces * 2, flush=True)

        return result

    def covariants(
            self,
            max_dimension,
            verbose=False,
            ignore_lower_dimension=False
    ):
        return self.cached_results.setdefault(
            ('covariants', max_dimension),
            self._covariants(max_dimension, verbose, ignore_lower_dimension)
        )

    @staticmethod
    def show_covariants(covariants):
        return "\n".join(
            "{irrep}: {operators}".format(
                irrep=irrep,
                operators=" + ".join(
                    "{count}{operator}{derivatives}".format(
                        count=str(count) + " " if count > 1 else "",
                        operator=operator,
                        derivatives="" if number_of_derivatives == 0 else (
                            " D" if number_of_derivatives == 1 else
                            f" D^{number_of_derivatives}"
                        )
                    )
                    for (operator, number_of_derivatives), count
                    in operator_counter.items()
                )
            )
            for irrep, operator_counter in covariants.items()
        )
