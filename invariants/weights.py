class Weight(object):
    def __init__(self, components):
        self.components = components

    def __str__(self):
        return "({})".format(
            " ".join(map(str, self))
        )

    __repr__ = __str__

    def __hash__(self):
        return hash(tuple(self))

    def __eq__(self, other):
        return self.components == other.components

    def __add__(self, other):
        return Weight([x + y for x, y in zip(self, other)])

    def __neg__(self):
        return Weight([-x for x in self])

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return Weight([other * x for x in self])

    __rmul__ = __mul__

    def __iter__(self):
        return iter(self.components)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return Weight(self.components[key])
        else:
            return self.components[key]

    def __len__(self):
        return len(self.components)

    @staticmethod
    def from_string(weight_str):
        return Weight(list(map(int, weight_str.split())))
