from invariants.algebras import Series, SimpleAlgebra, SemisimpleAlgebra
from invariants.representations import Irrep
from invariants.weights import Weight
from invariants.statistics import Statistics

from collections import Counter, defaultdict
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
            number_of_derivatives=0
    ):
        self.name = name
        self.lorentz_irrep = lorentz_irrep
        self.internal_irrep = internal_irrep
        self.charges = charges
        self.statistics = statistics
        self.dimension = dimension
        self.number_of_derivatives = number_of_derivatives

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return (
            self.name == other.name
            and self.number_of_derivatives == other.number_of_derivatives
            and self.irrep == other.irrep
        )

    def __str__(self):
        derivatives_str = (
            f"D^{self.number_of_derivatives} "
            if self.number_of_derivatives > 1
            else ("D" if self.number_of_derivatives == 1 else "")
        )
        return f"{derivatives_str}{self.name}"

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

    def differentiate(self, times):
        lorentz_irreps = itertools.chain.from_iterable(
            Irrep(_lorentz_algebra, Weight([n, n]))
            * self.lorentz_irrep
            for n in range(times % 2, times + 1, 2)
        )

        # !!
        lorentz_irreps = [
            Irrep(
                _lorentz_algebra,
                Weight([times, times]) + self.lorentz_irrep.highest_weight
            )
        ]

        return set(
            Field(
                name=self.name,
                lorentz_irrep=lorentz_irrep,
                internal_irrep=self.internal_irrep,
                charges=self.charges,
                statistics=self.statistics,
                dimension=self.dimension+times,
                number_of_derivatives=self.number_of_derivatives+times
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
    def irreps(self):
        power_irreps = (
            field.irrep.power(exponent, field.statistics)
            for field, exponent in self.content.items()
        )

        return Counter(functools.reduce(Irrep.product, power_irreps))

    @property
    def derivatives(self):
        results = set()
        for field in self.content:
            for new_field in field.differentiate(1):
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

    def internal_singlets(self, max_dimension):
        max_derivatives = int(max_dimension - self.dimension)

        return {
            number_of_derivatives: Counter(
                irrep
                for operator in self.differentiate(number_of_derivatives)
                for irrep in operator.irreps.elements()
                if irrep[2:].is_singlet
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

    def invariants(self, max_dimension):
        singlets = self.internal_singlets(max_dimension)
        print("SINGLETS:", singlets)

        for derivative_count in range(int(max_dimension - self.dimension)):
            print("count:", derivative_count)
            for initial_irrep, initial_count in singlets[derivative_count].items():
                descendants = Operator.descendants(
                    initial_irrep,
                    int(max_dimension - self.dimension - derivative_count),
                    int(derivative_count)
                )

                for number_of_derivatives, counter in descendants.items():
                    print(number_of_derivatives, counter)
                    singlets[number_of_derivatives] -= Counter({
                        irrep: count * initial_count
                        for irrep, count in counter.items()
                    })

            print("SINGLETS:", singlets)

        return {
            number_of_derivatives: sum(
                count for irrep, count in counter.items() if irrep.is_singlet
            )
            for number_of_derivatives, counter in singlets.items()
        }
            
                
                
            
        
# class Operator(object):
#     show_irrep = False

#     def __init__(self, ):
#         self.field_content = field_content

#     def __eq__(self, other):
#         return (
#             self.field_content == other.field_content
#             and self.lorentz_irrep == other.lorentz_irrep
#             and self.internal_irrep == other.internal_irrep
#         )

#     def __hash__(self):
#         return hash((
#             frozenset(self.field_content),
#             self.lorentz_irrep,
#             self.internal_irrep
#         ))

#     def __str__(self):
#         if len(self.field_content) == 0:
#             content = "1"
#         else:
#             content = " ".join(
#                 f"{field}^{exponent}" if exponent > 1 else f"{field}"
#                 for field, exponent in self.field_content.items()
#             )

#         if Operator.show_irrep:
#             return f"({content}: ({self.lorentz_irrep}, {self.internal_irrep}))"
#         else:
#             return content

#     __repr__ = __str__

#     @staticmethod
#     def _generate(field_content, irreps):
#         return Counter({
#             Operator(
#                 field_content,
#                 irrep.split(2)[0],
#                 irrep.split(2)[1]
#             ):
#             count
#             for irrep, count in irreps.items()
#         })

#     @property
#     def charges(self):
#         return [
#             sum(charges) for charges in zip(*[
#                 [
#                     charge * self.field_content[field]
#                     for charge in field.charges
#                 ]
#                 for field in self.field_content
#             ])
#         ]

#     @property
#     def dimension(self):
#         return sum(
#             field.dimension * self.field_content[field]
#             for field in self.field_content
#         )

#     @property
#     def is_invariant(self):
#         return (
#             all(charge == 0 for charge in self.charges)
#             and all(component == 0 for component in self.irrep.highest_weight)
#         )            


class EFT(object):
    def __init__(self, algebra, fields, use_eom=False):
        self.algebra = algebra
        self.fields = fields
        self.use_eom = use_eom

#     @staticmethod
#     def _combinations(fields, max_dimension):
#         if len(fields) == 0:
#             return [collections.Counter()]

#         max_exponent = math.floor(max_dimension / fields[0].dimension)

#         return [
#             collections.Counter({fields[0]: exponent}) + combination
#             for exponent in range(max_exponent + 1)
#             for combination in EFT._combinations(
#                     fields[1:],
#                     max_dimension - exponent * fields[0].dimension
#             )
#         ]

#     @property
#     def _singlet(self):
#         return Irrep.singlet(self.algebra)

#     def covariants(self, max_dimension):
#         operators = collections.Counter()

#         for field_content in EFT._combinations(self.fields, max_dimension):
#             powers = (
#                 field.irrep.power(exponent, field.statistics)
#                 for field, exponent in field_content.items()
#             )

#             irreps = functools.reduce(Irrep.product, powers, [self._singlet])

#             for irrep, count in collections.Counter(irreps).items():
#                 operators[Operator(field_content, irrep)] += count

#         return operators

#     @staticmethod
#     def _filter_invariants(operators):
#         return collections.Counter({
#             operator: count
#             for operator, count in operators.items()
#             if operator.is_invariant
#         })

#     def invariants(self, max_dimension):
#         return EFT._filter_invariants(
#             self.covariants(max_dimension)
#         )

# def differentiate_content(field_content):
#     out_contents = set()
#     for field in field_content:
#         out_content = field_content.clone()

#         out_content.add(field, -1)
#         out_content.add(field.derivative(), 1)

#         out_contents.add(out_content)

#     return out_contents


# def differentiate(algebra, operator):
#     vector_singlet = Irrep(algebra, Weight([1, 1] + [0] * algebra.rank))

#     return {
#         Operator(new_content, new_irrep)
#         for new_content in differentiate_content(operator.field_content)
#         for new_irrep in operator.irrep * vector_singlet
#     }


# def remove_one(dictionary):
#     for key in dictionary:
#         dictionary[key]


# def invariants_with_derivatives(algebra, field, max_dimension):
#     initial_operators = operators(algebra, field, max_dimension)
#     current_operators = initial_operators

#     all_invariants = filter_invariants(initial_operators)

#     for initial_operator in initial_operators:
#         current_operators = {initial_operator}

#         for number_of_derivatives in range(1, max_dimension):
#             # correct_irrep = (
#             #     initial_operator.irrep.highest_weight[:2]
#             #     == [number_of_derivatives, number_of_derivatives]
#             # )
#             # if not correct_irrep:
#             #     continue
            
#             current_operators = set(
#                 itertools.chain.from_iterable(
#                     differentiate(algebra, operator)
#                     for operator in current_operators
#                     if operator.dimension + 1 <= max_dimension
#                 )
#             )

#             new_invariants = {
#                 operator for operator in current_operators
#                 if operator.is_invariant()
#             }

#             for _ in range(number_of_derivatives):
#                 if new_invariants:
#                     new_invariants.pop()

#             all_invariants.update({
#                 operator.field_content: initial_operators[initial_operator]
#                 for operator in new_invariants
#             })

#     return all_invariants


def hilbert_series_str(operators):
    return " + ".join([
        f"{count} {op}" if count > 1 else str(op)
        for op, count in operators.items()
    ])
