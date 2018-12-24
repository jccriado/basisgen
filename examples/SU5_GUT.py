from basisgen import (
    irrep, algebra, boson, fermion,  scalar, L_spinor, R_spinor, Field, EFT
)

irrep_5 = irrep('SU5', '1 0 0 0')
irrep_10 = irrep('SU5', '0 1 0 0')
irrep_15 = irrep('SU5', '2 0 0 0')
irrep_24 = irrep('SU5', '1 0 0 1')

Phi = Field(
    name='Phi',
    lorentz_irrep=scalar,
    internal_irrep=irrep_24,
    statistics=boson,
    dimension=1
)

Hu = Field(
    name='Hu',
    lorentz_irrep=scalar,
    internal_irrep=irrep_5,
    statistics=boson,
    dimension=1
)

Hd = Field(
    name='Hd',
    lorentz_irrep=scalar,
    internal_irrep=irrep_5.conjugate,
    statistics=boson,
    dimension=1
)

psi = Field(
    name='psi',
    lorentz_irrep=L_spinor,
    internal_irrep=irrep_5.conjugate,
    statistics=fermion,
    dimension=1.5
)

Psi = Field(
    name='Psi',
    lorentz_irrep=L_spinor,
    internal_irrep=irrep_10,
    statistics=fermion,
    dimension=1.5
)

N = Field(
    name='N',
    lorentz_irrep=R_spinor,
    internal_irrep=irrep('SU5', '0 0 0 0'),
    statistics=fermion,
    dimension=1.5
)

fields = [Phi, Hu, Hd, psi, Psi, N]
SU5_GUT = EFT(algebra('SU5'), fields + [field.conjugate for field in fields])

print(SU5_GUT.invariants(max_dimension=4, verbose=True))
