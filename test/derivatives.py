from invariants.fields import Operator, EFT
from standard_model import *
from collections import Counter
import functools

# print(hilbert_series_str(
#     invariants_with_derivatives(SM_algebra, SM_fields, 6)
# ))

Operator.show_irrep = True
# print(list(
#     operator.derivative()
#     for operator in
#     Operator(Counter([phi, phic]), standard_model._singlet).derivative()
# ))

#print(Operator(Counter([phi, phic]), standard_model._singlet).derivative(singlet=standard_model._singlet))

# vector_singlet_irrep = Irrep(
#     SM_algebra,
#     Weight([
#         0 if index not in (0, 1) else 1
#         for index in range(5)
#     ])
# )

# irreps = functools.reduce(
#     Irrep.product,
#     [[vector_singlet_irrep]] * 3,
#     [phi.irrep]
# )

# print(irreps)

#print([f.irrep for f in phi.differentiate(3).elements()])

#print(EFT(SM_algebra, [phi, phic]).invariants_with_derivatives(6))
#result= standard_model.invariants_with_derivatives(6)
#print(result)
#print(len(result))

# print((phi*phic**2*dc*q).invariants(8))
invariants = SMEFT.invariants(6)
print("Number of invariants: {}".format(
    sum(
        count
        for counter in invariants.values()
        for count in counter.values()
    )
))
print(EFT.show_invariants(invariants))

