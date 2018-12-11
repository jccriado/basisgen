from invariants.weights import Weight
from invariants.statistics import Statistics
from invariants.multimap import MultivaluedMap

import collections
import copy
import itertools
import functools


class Representation(object):
    def __init__(self, weights):
        self.weights = weights

    def __str__(self):
        return str(self.weights)

    __repr__ = __str__

    def __iter__(self):
        return iter(self.weights)

    def __mul__(self, other):
        return Representation(collections.Counter(
            first_weight + second_weight
            for first_weight in self.weights.elements()
            for second_weight in other.weights.elements()
        ))

    @functools.lru_cache(maxsize=None)
    def highest_weight(self, algebra):
        return max(self.weights, key=algebra.height)

    @functools.lru_cache(maxsize=None)
    def decompose(self, algebra):
        remaining_weights = copy.copy(self.weights)
        irreps = []

        while remaining_weights:
            highest_weight = Representation(
                remaining_weights
            ).highest_weight(algebra)

            current_irrep = Irrep(algebra, highest_weight)
            irreps.append(current_irrep)

            remaining_weights -= current_irrep.representation.weights

        return irreps


class Irrep(object):
    def __init__(self, algebra, highest_weight):
        self.algebra = algebra
        self.highest_weight = highest_weight
        if isinstance(self.algebra, list):
            raise Exception(str(algebra))

    def __str__(self):
        return "Irrep({algebra}, {highest_weight})".format(
            algebra=self.algebra,
            highest_weight=self.highest_weight
        )

    __repr__ = __str__

    def __add__(self, other):
        concatenated_weight = Weight(list(
            itertools.chain(self.highest_weight, other.highest_weight)
        ))
        return Irrep(self.algebra + other.algebra, concatenated_weight)

    def __hash__(self):
        return hash((self.algebra, tuple(self.highest_weight)))

    def __eq__(self, other):
        return (
            self.algebra == other.algebra
            and tuple(self.highest_weight) == tuple(other.highest_weight)
        )

    @functools.lru_cache(maxsize=None)
    def __mul__(self, other):
        if isinstance(other, Irrep):
            product_representation = self.representation * other.representation
            return product_representation.decompose(self.algebra)
        else:
            return Irrep.product([self], other)

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

    @staticmethod
    def product(first_irrep_list, second_irrep_list):
        return list(
            itertools.chain.from_iterable([
                first_irrep * second_irrep
                for first_irrep in first_irrep_list
                for second_irrep in second_irrep_list
            ])
        )

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
            for alpha in self.algebra.positive_roots
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
    def _semisimple_representation(self):
        algebras = self.algebra.simple_algebras
        split_highest_weight = self.algebra.split_weight(self.highest_weight)

        split_weights = (
            Irrep(algebra, weight).representation.weights.elements()
            for algebra, weight in zip(algebras, split_highest_weight)
        )

        return Representation(collections.Counter([
            Weight(list(itertools.chain.from_iterable(weights)))
            for weights in itertools.product(*split_weights)
        ]))

    @property
    @functools.lru_cache(maxsize=None)
    def representation(self):
        if isinstance(self.algebra, collections.Iterable):
            return self._semisimple_representation
        else:
            return Representation(self.weights_with_multiplicities)

    @functools.lru_cache(maxsize=None)
    def power(self, exponent, statistics):
        combinations_function = {
            Statistics.BOSON: itertools.combinations_with_replacement,
            Statistics.FERMION: itertools.combinations
        }[statistics]

        weights = self.representation.weights.elements()

        power_weights = collections.Counter([
            sum(combination, Irrep.singlet(self.algebra).highest_weight)
            for combination in combinations_function(weights, exponent)
        ])

        return Representation(power_weights).decompose(self.algebra)


    # def _weight_multiplicity(self, level, weight, previous_multiplicities):
    #     if level == 0:
    #         return 1

    #     delta = self.algebra.sum_of_positive_roots

    #     weight = np.array(weight.components)
    #     positive_roots = np.array(self.algebra.positive_roots)
    #     k_values = np.array(range(1, level + 1))
    #     metric = np.array(self.algebra.metric)

    #     weights_plus_k_alpha = (
    #         weight[np.newaxis, np.newaxis, ...]
    #         + positive_roots[np.newaxis, ...]
    #         * k_values[..., np.newaxis, np.newaxis]
    #     )

    #     scalar_products = np.einsum(
    #         'ijk,kl,jl->ij',
    #         weights_plus_k_alpha,
    #         metric,
    #         positive_roots,
    #         optimize=True
    #     )

    #     numerator = 2 * sum(
    #         previous_multiplicities.get(Weight(list(weight_plus_k_alpha)), 0)
    #         * scalar_product
    #         for scalar_product, weight_plus_k_alpha in zip(
    #                 scalar_products.reshape((-1,)),
    #                 weights_plus_k_alpha.reshape((-1, len(weight)))
    #         )
    #     )

    #     denominator = (
    #         self.algebra.norm_squared(self.highest_weight + delta)
    #         - self.algebra.norm_squared(Weight(list(weight)) + delta)
    #     )

    #     return int(round(numerator / denominator))
