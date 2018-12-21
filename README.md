[![Build Status](https://travis-ci.com/jccriado/invariants.svg?token=FCho83xJ9cZySjbkvBWS&branch=master)](https://travis-ci.com/jccriado/invariants)

This is a package for calculating bases of operators in effective field
theories. It can also be used to compute independent sets of
gauge-covariant operators. It includes functionality for doing calculations with
semisimple Lie algebra representations.

To obtain a basis of operators, the following input is needed:
- A semisimple Lie algebra. This corresponds to the symmetry group of the
  field theory.
- A set of fields, with information about their irreducible representation
  (irrep) under the symmetry group, their charges under an arbitrary number of
  _U(1)_ factors, their dimension, statistics, etc. 

Fixing a maximum dimension, a set of invariant operators constructed using the
fields and their derivatives can be computed. Integration by parts redundancy
is automatically taken into account. Optionally, equations of motion can be
imposed to reduce the number of invariants.

## Installation

Just do `pip install invariants`. Python 3.6+ is needed.

## Usage

### Working with semisimple Lie algebra representations

As an example, we will compute the weight system for the representation of
A<sub>2</sub> = su(3) with highest weight `(1 1)` (an octet):

~~~
>>> from invariants.shortcuts import irrep
>>> irrep('A2', '1 1').weights_view()
    (1 1)
(2 -1) (-1 2)
 (0 0) (0 0)
(1 -2) (-2 1)
   (-1 -1)
~~~

The decomposition of the tensor product of `(1 0)` and `(0 1)` (a triplet
and an anti-triplet) is done directly as:

~~~
>>> irrep('A2', '1 0') * irrep('A2', '0 1')
[(1 1)] + [(0 0)]
~~~

This also works for semisimple algebras:

~~~
>>> irrep('B2 + G2', '0 2 1 0 0 0').weights_view()
     (0 2 1 0)
     (0 2 -1 3)
      (0 2 0 1)
     (0 2 1 -1)
(0 2 -1 2) (0 2 2 -3)
 (0 2 0 0) (0 2 0 0)
(0 2 1 -2) (0 2 -2 3)
     (0 2 -1 1)
     (0 2 0 -1)
     (0 2 1 -3)
     (0 2 -1 0)
~~~

The group notation (using `SO`, `SU` and `Sp` for simple algebras and `x` for
their direct sum) can be used. For example an _SU(3)_ triplet, _SU(2)_ doublet
would be:

~~~
>>> irrep('SU3 x SU2', '1 0 1').weights_view()
     (1 0 1)
(-1 1 1) (1 0 -1)
(0 -1 1) (-1 1 -1)
    (0 -1 -1)
~~~

### EFT operators

#### Simple example

Consider a simple EFT: a scalar _SU(2)_ doublet with hypercharge 1/2. To
compute all the independent invariants built with this field and its derivatives
(using integration by parts and equations of motion) we can write the script:

~~~
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
~~~

The output is

~~~
phi phi*: 1
(phi)^2 (phi*)^2: 1
(phi)^2 (phi*)^2 D^2: 2
(phi)^2 (phi*)^2 D^4: 3
(phi)^3 (phi*)^3: 1
(phi)^3 (phi*)^3 D^2: 2
(phi)^4 (phi*)^4: 1
Total: 11
~~~
	
Each line gives the number of independent invariant operators with the
specified field content (and number covariant derivatives).

#### The Standard Model EFT

The SMEFT is defined in `invariants.smeft`. An example of use of this module can
be found in `examples/standard_model.py`. Runnning it as

~~~
python standard_model.py --dimension 8
~~~

Gives 1123 invariants (in ~ 40 seconds in a 2,6 GHz Intel Core i5). This agrees
with arXiv:1512.03433.
Ignoring lower dimensional operators in each case:
* Dimension 9: 560 (3 minutes)
* Dimension 10: 15456 (15 minutes)
* 


