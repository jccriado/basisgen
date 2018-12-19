from invariants.eft import Field, EFT
from invariants.statistics import Statistics

from invariants.shortcuts import (
    SU2_algebra, SU2_singlet, SU2_doublet, SU2_triplet,
    SU3_algebra, SU3_singlet, SU3_triplet, SU3_anti_triplet, SU3_octet,
    scalar, L_spinor, R_spinor, L_tensor, R_tensor
)

sm_gauge_algebra = SU3_algebra + SU2_algebra


phi = Field(
    name='phi',
    lorentz_irrep=scalar,
    internal_irrep=SU3_singlet+SU2_doublet,
    charges=[1/2],
    statistics=Statistics.BOSON,
    dimension=1
)

phic = Field(
    name='phic',
    lorentz_irrep=scalar,
    internal_irrep=SU3_singlet+SU2_doublet,
    charges=[-1/2],
    statistics=Statistics.BOSON,
    dimension=1
)

bL = Field(
    name='bL',
    lorentz_irrep=L_tensor,
    internal_irrep=SU3_singlet+SU2_singlet,
    charges=[0],
    statistics=Statistics.BOSON,
    dimension=2
)

bR = Field(
    name='bR',
    lorentz_irrep=R_tensor,
    internal_irrep=SU3_singlet+SU2_singlet,
    charges=[0],
    statistics=Statistics.BOSON,
    dimension=2
)

wL = Field(
    name='wL',
    lorentz_irrep=L_tensor,
    internal_irrep=SU3_singlet+SU2_triplet,
    charges=[0],
    statistics=Statistics.BOSON,
    dimension=2
)

wR = Field(
    name='wR',
    lorentz_irrep=R_tensor,
    internal_irrep=SU3_singlet+SU2_triplet,
    charges=[0],
    statistics=Statistics.BOSON,
    dimension=2
)

gL = Field(
    name='gL',
    lorentz_irrep=L_tensor,
    internal_irrep=SU3_octet+SU2_singlet,
    charges=[0],
    statistics=Statistics.BOSON,
    dimension=2
)

gR = Field(
    name='gR',
    lorentz_irrep=R_tensor,
    internal_irrep=SU3_octet+SU2_singlet,
    charges=[0],
    statistics=Statistics.BOSON,
    dimension=2
)


def q(number_of_flavors=1):
    return Field(
        name='q',
        lorentz_irrep=L_spinor,
        internal_irrep=SU3_triplet+SU2_doublet,
        charges=[1/6],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def qc(number_of_flavors=1):
    return Field(
        name='qc',
        lorentz_irrep=R_spinor,
        internal_irrep=SU3_anti_triplet+SU2_doublet,
        charges=[-1/6],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def u(number_of_flavors=1):
    return Field(
        name='u',
        lorentz_irrep=R_spinor,
        internal_irrep=SU3_triplet+SU2_singlet,
        charges=[2/3],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def uc(number_of_flavors=1):
    return Field(
        name='uc',
        lorentz_irrep=L_spinor,
        internal_irrep=SU3_anti_triplet+SU2_singlet,
        charges=[-2/3],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def d(number_of_flavors=1):
    return Field(
        name='d',
        lorentz_irrep=R_spinor,
        internal_irrep=SU3_triplet+SU2_singlet,
        charges=[-1/3],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def dc(number_of_flavors=1):
    return Field(
        name='dc',
        lorentz_irrep=L_spinor,
        internal_irrep=SU3_anti_triplet+SU2_singlet,
        charges=[1/3],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def L(number_of_flavors=1):
    return Field(
        name='l',
        lorentz_irrep=L_spinor,
        internal_irrep=SU3_singlet+SU2_doublet,
        charges=[-1/2],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def Lc(number_of_flavors=1):
    return Field(
        name='lc',
        lorentz_irrep=R_spinor,
        internal_irrep=SU3_singlet+SU2_doublet,
        charges=[1/2],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def e(number_of_flavors=1):
    return Field(
        name='e',
        lorentz_irrep=R_spinor,
        internal_irrep=SU3_singlet+SU2_singlet,
        charges=[-1],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


def ec(number_of_flavors=1):
    return Field(
        name='ec',
        lorentz_irrep=L_spinor,
        internal_irrep=SU3_singlet+SU2_singlet,
        charges=[1],
        statistics=Statistics.FERMION,
        dimension=1.5,
        number_of_flavors=number_of_flavors
    )


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
