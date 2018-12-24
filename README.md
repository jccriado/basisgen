[![Build Status](https://travis-ci.com/jccriado/basisgen.svg?token=FCho83xJ9cZySjbkvBWS&branch=master)](https://travis-ci.com/jccriado/basisgen)

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

Just do `pip install basisgen`. Python 3.5+ is needed.

## Usage examples

### Working with semisimple Lie algebra representations

As an example, we will compute the weight system for the representation of
A<sub>2</sub> = su(3) with highest weight `(1 1)` (an octet):

``` python-console
>>> from basisgen import irrep
>>> irrep('A2', '1 1').weights_view()
    (1 1)
(2 -1) (-1 2)
 (0 0) (0 0)
(1 -2) (-2 1)
   (-1 -1)
```

The decomposition of the tensor product of `(1 0)` and `(0 1)` (a triplet
and an anti-triplet) is done directly as:

``` python-console
>>> irrep('A2', '1 0') * irrep('A2', '0 1')
[(1 1)] + [(0 0)]
```

This also works for semisimple algebras:

``` python-console
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
```

The group notation (using `SO`, `SU` and `Sp` for simple algebras and `x` for
their direct sum) can be used. For example an _SU(3)_ triplet, _SU(2)_ doublet
would be:

``` python-console
>>> irrep('SU3 x SU2', '1 0 1').weights_view()
     (1 0 1)
(-1 1 1) (1 0 -1)
(0 -1 1) (-1 1 -1)
    (0 -1 -1)
```

### EFT operators

#### Simple example

Consider a simple EFT: a scalar _SU(2)_ doublet with hypercharge 1/2. To
compute all the independent invariants built with this field and its derivatives
(using integration by parts and equations of motion) we can write the script:

```python
from basisgen import algebra, irrep, scalar, Field, EFT

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
```

This code can be also found in `examples/simple.py`. Running it gives:

``` shell
$ python simple.py
phi phi*: 1
(phi)^2 (phi*)^2: 1
(phi)^2 (phi*)^2 D^2: 2
(phi)^2 (phi*)^2 D^4: 3
(phi)^3 (phi*)^3: 1
(phi)^3 (phi*)^3 D^2: 2
(phi)^4 (phi*)^4: 1
Total: 11
```
	
Each line gives the number of independent invariant operators with the
specified field content (and number covariant derivatives).

#### The Standard Model EFT

The SMEFT is defined in `basisgen.smeft`. See the code there for a more complex
example. The script `examples/standard_model.py` makes use of this module. To
obtain the dimension 8 operators in the SMEFT (with one generation), do:

``` shell
$ python standard_model.py --dimension 8
```

This gives 993 invariants (in ~ 40 seconds in a 2,6 GHz Intel Core i5). This
agrees with arXiv:1512.03433. For dimension 9, it gives 560 operators (in 3
minutes), and for dimension 10, 15456 (in 15 minutes).

Do `python standard_model.py --help` to see more options included in the script.
For example, the dimension-6 SMEFT with 3 generations of fermions is obtained by:

``` shell
$ python standard_model.py --dimension 6 --number_of_flavors 3
```


#### SU(5) GUT example

As another example, consider the Georgi-Glashow model of grand unification.  The
internal symmetry group is _SU(5)_. In `examples/SU5_GUT.py` the code that
defines this theory can be found. It outputs 83 possible field contents for
operators, each of them with corresponding the number of independent operators
of that form.

## Summary of basic objects

### Functions

- `algebra`. Creates a (semi)simple Lie algebra (`algebras.SimpleAlgebra` or
  `algebras.SemisimpleAlgebra`) from one string argument. Examples of valid
  arguments are: `'A3'`, `'C12'`, `'F4'`, `'SU3'`, `'B2 + E7'` or `'SU5 x SO6 x
  Sp10'`.

- `irrep`. Creates an irreducible representation (`representations.Irrep`) from
  2 string arguments: the first represents the algebra and the second the
  highest weight. Example: `irrep('SU4 x Sp7', '1 0 1 0 2 1')`. The weight
  system of the an irrep can be shown by calling its `weights_view()`
  method. Irreps can be multiplied to get the decomposition of their tensor
  product.

### Classes

- `Field`. Useful attribute: `conjugate`. Constructor arguments:

| Name                | Description                                         | Default |
|---------------------|:----------------------------------------------------|---------|
| `name`              | String Identifier                                   |         |
| `lorentz_irrep`     | Lorentz group irrep                                 |         |
| `internal_irrep`    | Irrep of the internal symmetry group of the theory  |         |
| `charges`           | Charges under an arbitrary number of _U(1)_ factors | `[]`    |
| `statistics`        | Either `boson` or `fermion`                         | `boson` |
| `dimension`         | Energy dimensions of the field                      | `1`     |
| `number_of_flavors` | Number of different copies of the same field        | `1`     |


- `EFT`. Constructor arguments: 
  - `internal_algebra`, the semisimple Lie algebra of the internal symmetry group
  - `fields`, a list of `Field` objects representing the field content of the theory.
  
  Methods:
  - `invariants(max_dimension, verbose=False)`. Return a basis of operators,
  encapsulated in an `EFT.Invariants` object.
  - `covariants(max_dimension, verbose=False)`. Returns a collection of all
  operators with all possible irreps, in the form of a `EFT.Covariants`
  instance.
  
  Classes:
  - `EFT.Invariants`. Can be directly printed (implements `__str__`).  Has a
  method `count` to calculate the total number of operators in the basis. Has a
  method `show_by_classes`, which returns a simplified string representation of
  the basis, provided a dictionary whose keys are the fields and values are
  strings representing classes of fields.
  - `EFT.Covariants`. Its only purpose is to hold the information until it is
  printed (implements `__str__`).

### Variables

- Predefined irreps of the Lorentz group: `scalar`, `L_spinor`, `R_spinor`,
  `vector`, `L_tensor`, `R_tensor`.
- Predefined statistics: `boson`, `fermion`.
