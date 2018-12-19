from invariants.eft import Field, EFT
from invariants.shortcuts import irrep, algebra, scalar

phi = Field(
    name='phi',
    lorentz_irrep=scalar,
    internal_irrep=irrep('SU2', '1'),
    charges=[1/2]
)
phic = Field(
    name='phi*',
    lorentz_irrep=scalar,
    internal_irrep=irrep('SU2', '1'),
    charges=[-1/2]
)

my_eft = EFT(algebra('SU2'), [phi, phic])

invariants = my_eft.invariants(max_dimension=8)

print(invariants)
print("Total:", invariants.count())
