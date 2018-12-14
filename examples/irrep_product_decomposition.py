from invariants.SU3 import SU3_irrep, SU3_show


if __name__ == '__main__':
    result = SU3_irrep('3') * SU3_irrep('3*')
    print(f"3 x 3* = " + SU3_show(result))

    result = SU3_irrep('3') * SU3_irrep('6')
    print(f"3 x 6 = " + SU3_show(result))

    result = SU3_irrep('3*') * SU3_irrep('8') * SU3_irrep('3')
    print(f"3* x 8 x 3 = " + SU3_show(result))
