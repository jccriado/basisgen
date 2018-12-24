from invariants import algebra, irrep, scalar, Field, EFT

phi = Field(
    name='phi',
    lorentz_irrep=scalar,
    internal_irrep=irrep('SU2', '1'),
    charges=[1/2]
)

my_eft = EFT(algebra('SU2'), [phi, phi.conjugate])

invariants = my_eft.invariants(max_dimension=8)

print(invariants)
print("Total:", invariants.count())
