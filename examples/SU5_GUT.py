from basisgen import irrep, algebra, boson, scalar, Field, EFT

irrep_5 = irrep('SU5', '1 0 0 0')
irrep_24 = irrep('SU5', '1 0 0 1')

Phi = Field(
    name='Phi',
    lorentz_irrep=scalar,
    internal_irrep=irrep_24,
    statistics=boson,
    dimension=1
)

phi = Field(
    name='phi',
    lorentz_irrep=scalar,
    internal_irrep=irrep_5,
    statistics=boson,
    dimension=1
)

SU5_GUT_scalar_sector = EFT(algebra('SU5'), [Phi, phi, phi.conjugate])

print(SU5_GUT_scalar_sector.invariants(max_dimension=4, verbose=True))
