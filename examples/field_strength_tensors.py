from basisgen import irrep, algebra, Field, EFT

BL, BR = Field.strength_tensors(
    name='B',
    internal_irrep=irrep('SU2', '0'),
    charges=[0],
)

my_eft = EFT(algebra('SU2'), [BL, BR])

invariants = my_eft.invariants(max_dimension=8, use_eom=False)

print(invariants)
