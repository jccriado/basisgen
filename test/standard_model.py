from invariants.weights import Weight
from invariants.algebras import SimpleAlgebra
from invariants.statistics import Statistics
from invariants.representations import Irrep
from invariants.fields import Field, EFT

SU2 = SimpleAlgebra.from_string('A1')
SU3 = SimpleAlgebra.from_string('A2')
Lorentz = SU2 + SU2

SU2_singlet = Irrep(SU2, Weight([0]))
SU2_doublet = Irrep(SU2, Weight([1]))
SU2_triplet = Irrep(SU2, Weight([2]))

SU3_singlet = Irrep(SU3, Weight([0, 0]))
SU3_triplet = Irrep(SU3, Weight([1, 0]))
SU3_antitriplet = Irrep(SU3, Weight([0, 1]))
SU3_octet = Irrep(SU3, Weight([1, 1]))

scalar = Irrep(Lorentz, Weight([0, 0]))
L_spinor = Irrep(Lorentz, Weight([1, 0]))
R_spinor = Irrep(Lorentz, Weight([0, 1]))
vector = Irrep(Lorentz, Weight([1, 1]))
L_tensor = Irrep(Lorentz, Weight([2, 0]))
R_tensor = Irrep(Lorentz, Weight([0, 2]))

SM_algebra = Lorentz + SU3 + SU2

q = Field(
    'q', L_spinor, SU3_triplet + SU2_doublet, [1/6], Statistics.FERMION, 1.5
)
qc = Field(
    'qc', R_spinor, SU3_antitriplet + SU2_doublet, [-1/6], Statistics.FERMION,
    1.5
)

u = Field(
    'u', R_spinor, SU3_triplet + SU2_singlet, [2/3], Statistics.FERMION, 1.5
)
uc = Field(
    'uc', L_spinor, SU3_antitriplet + SU2_singlet, [-2/3], Statistics.FERMION,
    1.5
)

d = Field(
    'd', R_spinor, SU3_triplet + SU2_singlet, [-1/3], Statistics.FERMION, 1.5
)
dc = Field(
    'dc', L_spinor, SU3_antitriplet + SU2_singlet, [1/3], Statistics.FERMION,
    1.5
)

l = Field(
    'l', L_spinor, SU3_singlet + SU2_doublet, [-1/2], Statistics.FERMION, 1.5
)
lc = Field(
    'lc', R_spinor, SU3_singlet + SU2_doublet, [1/2], Statistics.FERMION, 1.5
)

e = Field(
    'e', R_spinor, SU3_singlet + SU2_singlet, [-1], Statistics.FERMION, 1.5
)
ec = Field(
    'ec', L_spinor, SU3_singlet + SU2_singlet, [1], Statistics.FERMION, 1.5
)

phi = Field(
    'phi', scalar, SU3_singlet + SU2_doublet, [1/2], Statistics.BOSON, 1
)
phic = Field(
    'phic', scalar, SU3_singlet + SU2_doublet, [-1/2], Statistics.BOSON, 1
)

bL = Field(
    'bL', L_tensor, SU3_singlet + SU2_singlet, [0], Statistics.BOSON, 2
)
bR = Field(
    'bR', R_tensor, SU3_singlet + SU2_singlet, [0], Statistics.BOSON, 2
)

wL = Field(
    'wL', L_tensor, SU3_singlet + SU2_triplet, [0], Statistics.BOSON, 2
)
wR = Field(
    'wR', R_tensor, SU3_singlet + SU2_triplet, [0], Statistics.BOSON, 2
)

gL = Field(
    'gL', L_tensor, SU3_octet + SU2_singlet, [0], Statistics.BOSON, 2
)
gR = Field(
    'gR', R_tensor, SU3_octet + SU2_singlet, [0], Statistics.BOSON, 2
)

SM_fields = [
    q, qc,
    u, uc, d, dc,
    l, lc,
    e, ec,
    phi, phic,
    bL, bR, wL, wR, gL, gR
]

SMEFT = EFT(SM_algebra, SM_fields)

if __name__ == "__main__":
    print(SMEFT.invariants(6))
