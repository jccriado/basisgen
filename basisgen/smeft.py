from basisgen import (
    irrep, algebra, scalar, L_spinor, R_spinor, L_tensor,
    boson, fermion, Field, EFT
)

sm_gauge_algebra = algebra('SU3 x SU2')


def sm_irrep(highest_weight_str):
    return irrep('SU3 x SU2', highest_weight_str)


phi = Field(
    name='phi',
    lorentz_irrep=scalar,
    internal_irrep=sm_irrep('0 0 1'),
    charges=[1/2],
    statistics=boson,
    dimension=1
)
phic = phi.conjugate

bL = Field(
    name='bL',
    lorentz_irrep=L_tensor,
    internal_irrep=sm_irrep('0 0 0'),
    charges=[0],
    statistics=boson,
    dimension=2
)
bR = bL.conjugate


wL = Field(
    name='wL',
    lorentz_irrep=L_tensor,
    internal_irrep=sm_irrep('0 0 2'),
    charges=[0],
    statistics=boson,
    dimension=2
)
wR = wL.conjugate

gL = Field(
    name='gL',
    lorentz_irrep=L_tensor,
    internal_irrep=sm_irrep('1 1 0'),
    charges=[0],
    statistics=boson,
    dimension=2
)
gR = gL.conjugate


def q(number_of_flavors=1):
    return Field(
        name='q',
        lorentz_irrep=L_spinor,
        internal_irrep=sm_irrep('1 0 1'),
        charges=[1/6],
        statistics=fermion,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def qc(number_of_flavors=1):
    return q(number_of_flavors).conjugate


def u(number_of_flavors=1):
    return Field(
        name='u',
        lorentz_irrep=R_spinor,
        internal_irrep=sm_irrep('1 0 0'),
        charges=[2/3],
        statistics=fermion,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def uc(number_of_flavors=1):
    return u(number_of_flavors).conjugate


def d(number_of_flavors=1):
    return Field(
        name='d',
        lorentz_irrep=R_spinor,
        internal_irrep=sm_irrep('1 0 0'),
        charges=[-1/3],
        statistics=fermion,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def dc(number_of_flavors=1):
    return d(number_of_flavors).conjugate


def L(number_of_flavors=1):
    return Field(
        name='l',
        lorentz_irrep=L_spinor,
        internal_irrep=sm_irrep('0 0 1'),
        charges=[-1/2],
        statistics=fermion,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def Lc(number_of_flavors=1):
    return L(number_of_flavors).conjugate


def e(number_of_flavors=1):
    return Field(
        name='e',
        lorentz_irrep=R_spinor,
        internal_irrep=sm_irrep('0 0 0'),
        charges=[-1],
        statistics=fermion,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def ec(number_of_flavors=1):
    return e(number_of_flavors).conjugate


def sm_fermions(number_of_flavors):
    return [
        fermion(number_of_flavors)
        for fermion in [q, qc, u, uc, d, dc, L, Lc, e, ec]
    ]


sm_field_strengths = [bL, bR, wL, wR, gL, gR]

sm_scalars = [phi, phic]


def sm_field_classes(number_of_flavors=1):
    out = {}
    out.update({field: 'phi' for field in sm_scalars})
    out.update({field: 'F' for field in sm_field_strengths})
    out.update({field: 'psi' for field in sm_fermions(1)})

    return out


def smeft(number_of_flavors=1):
    return EFT(
        sm_gauge_algebra,
        sm_scalars + sm_field_strengths + sm_fermions(number_of_flavors)
    )
