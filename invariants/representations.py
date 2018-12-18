from invariants.algebras import SemisimpleAlgebra
from invariants.weights import Weight
from invariants.statistics import Statistics
from invariants.containers import MultivaluedMap, OrderedCounter

import collections
import itertools
import functools


class WeightSystem(object):
    def __init__(self, weights):
        self.weights = collections.Counter(weights)

    def __str__(self):
        return str(self.weights.elements())

    __repr__ = __str__

    def __iter__(self):
        return iter(self.weights.items())

    def __add__(self, other):
        return WeightSystem(self.weights + other.weights)

    def __mul__(self, other):
        return WeightSystem(
            first_weight + second_weight
            for first_weight in self.weights.elements()
            for second_weight in other.weights.elements()
        )

    @functools.lru_cache(maxsize=None)
    def highest_weight(self, algebra):
        return max(self.weights, key=algebra.height)

    def _sorted_weights(self, algebra):
        def _first_height(pair):
            return algebra.height(pair[0])

        return OrderedCounter.sort(self, _first_height)

    @functools.lru_cache(maxsize=None)
    def decompose(self, algebra):
        remaining_weights = self._sorted_weights(algebra)
        irreps = IrrepCounter()

        while remaining_weights:
            highest_weight, multiplicity = remaining_weights.popitem()

            current_irrep = Irrep(algebra, highest_weight)
            irreps[current_irrep] += multiplicity

            current_rep = current_irrep.weight_system
            remaining_weights[highest_weight] = multiplicity
            remaining_weights -= collections.Counter({
                weight: count * multiplicity
                for weight, count in current_rep
            })

        return irreps


class Irrep(object):
    def __init__(self, algebra, highest_weight):
        self.algebra = algebra
        self.highest_weight = highest_weight

    def __str__(self):
        return "[{weight}]".format(weight=self.highest_weight)

    def __repr__(self):
        return "Irrep({algebra}, {weight})".format(
            algebra=self.algebra,
            weight=self.highest_weight
        )

    def __add__(self, other):
        return Irrep(
            self.algebra + other.algebra,
            self.highest_weight.concat(other.highest_weight)
        )

    def __hash__(self):
        return hash((self.algebra, self.highest_weight))

    def __eq__(self, other):
        return (
            self.algebra == other.algebra
            and self.highest_weight == other.highest_weight
        )

    @staticmethod
    @functools.lru_cache(maxsize=None)
    def positive_roots(algebra):
        all_roots = Irrep(algebra, algebra.highest_root).weights_by_level

        return [
            root
            for level in range(algebra.level_of_simple_roots + 1)
            for root in all_roots[level]
        ]

    def split(self):
        return map(
            Irrep,
            self.algebra.simple_algebras,
            self.algebra.split_weight(self.highest_weight)
        )

    @functools.lru_cache(maxsize=None)
    def _mul_semisimple_irreps(self, other):
        out = IrrepCounter()

        zipped_with_product = (
            (first_irrep * second_irrep).items()
            for first_irrep, second_irrep
            in zip(self.split(), other.split())
        )

        for combination in itertools.product(*zipped_with_product):
            irrep, count = combination[0]
            for inner_irrep, inner_count in combination[1:]:
                irrep += inner_irrep
                count *= inner_count

            out += IrrepCounter({irrep: count})

        return out

    @functools.lru_cache(maxsize=None)
    def _mul_simple_irreps(self, other):
        product_weight_system = self.weight_system * other.weight_system
        return product_weight_system.decompose(self.algebra)

    def __mul__(self, other):
        if isinstance(other, Irrep):
            if isinstance(self.algebra, SemisimpleAlgebra):
                return self._mul_semisimple_irreps(other)
            else:
                return self._mul_simple_irreps(other)

        else:
            return Irrep.product(IrrepCounter([self]), other)

    __rmul__ = __mul__

    def __getitem__(self, index):
        return Irrep(self.algebra[index], self.highest_weight[index])

    @property
    def is_singlet(self):
        return all(component == 0 for component in self.highest_weight)

    @staticmethod
    @functools.lru_cache(maxsize=None)
    def singlet(algebra):
        return Irrep(algebra, Weight([0] * algebra.rank))

    def _child_weights(self, weight, level):
        return MultivaluedMap.from_pairs(
            (level + k, weight - k*root)
            for root, component in zip(self.algebra.simple_roots, weight)
            for k in range(1, component + 1)
        )

    @property
    def weights_by_level(self):
        current_weights = self._child_weights(self.highest_weight, 0)
        all_weights = MultivaluedMap({0: {self.highest_weight}})

        while current_weights:
            previous_weights = current_weights
            current_weights = MultivaluedMap()

            for level in previous_weights:
                for weight in previous_weights[level]:
                    current_weights.update(self._child_weights(weight, level))

            all_weights.update(previous_weights)

        return all_weights

    def _weight_multiplicity(self, level, weight, previous_multiplicities):
        if level == 0:
            return 1

        delta = self.algebra.sum_of_positive_roots

        numerator = 2 * sum(
            previous_multiplicities.get(weight + k*alpha, 0)
            * self.algebra.scalar_product(weight + k*alpha, alpha)
            for k in range(1, level + 1)
            for alpha in Irrep.positive_roots(self.algebra)
        )

        denominator = (
            self.algebra.norm_squared(self.highest_weight + delta)
            - self.algebra.norm_squared(weight + delta)
        )

        return round(numerator / denominator)

    @property
    def weights_with_multiplicities(self):
        multiplicities = collections.Counter()

        for level, weights in self.weights_by_level.items():
            for weight in weights:
                multiplicities[weight] = self._weight_multiplicity(
                    level,
                    weight,
                    multiplicities
                )

        return multiplicities

    @property
    def _semisimple_weight_system(self):
        algebras = self.algebra.simple_algebras
        split_highest_weight = self.algebra.split_weight(self.highest_weight)

        split_weights = (
            Irrep(algebra, weight).weight_system.weights.elements()
            for algebra, weight in zip(algebras, split_highest_weight)
        )

        return WeightSystem(
            Weight(list(itertools.chain.from_iterable(weights)))
            for weights in itertools.product(*split_weights)
        )

    @property
    def weight_system(self):
        if isinstance(self.algebra, SemisimpleAlgebra):
            return self._semisimple_weight_system
        else:
            return WeightSystem(self.weights_with_multiplicities)

    @functools.lru_cache(maxsize=None)
    def power(self, exponent, statistics):
        combinations_function = {
            Statistics.BOSON: itertools.combinations_with_replacement,
            Statistics.FERMION: itertools.combinations
        }[statistics]

        weights = self.weight_system.weights.elements()

        power_weights = (
            sum(combination, Irrep.singlet(self.algebra).highest_weight)
            for combination in combinations_function(weights, exponent)
        )

        return WeightSystem(power_weights).decompose(self.algebra)


class IrrepCounter(collections.Counter):
    def __str__(self):
        def term_str(irrep, counter):
            return "{counter}{irrep}".format(
                counter="" if counter == 1 else str(counter) + " ",
                irrep=str(irrep)
            )

        return " + ".join(
            term_str(irrep, counter) for irrep, counter in self.items()
        )

    def __add__(self, other):
        return IrrepCounter(super().__add__(other))

    def __mul__(self, other):
        if not isinstance(other, collections.Counter):
            other = IrrepCounter([other])

        return sum(
            (
                IrrepCounter({
                    irrep: count * first_count * second_count
                    for irrep, count in (first_irrep * second_irrep).items()
                })
                for first_irrep, first_count in self.items()
                for second_irrep, second_count in other.items()
            ),
            IrrepCounter()
        )
