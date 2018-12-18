class Weight(object):
    def __init__(self, components):
        self.components = tuple(components)

    def __str__(self):
        return "({})".format(
            " ".join(map(str, self))
        )

    __repr__ = __str__

    def __hash__(self):
        return hash(self.components)

    def __eq__(self, other):
        return self.components == other.components

    def __add__(self, other):
        return Weight(x + y for x, y in zip(self, other))

    def __neg__(self):
        return Weight(-x for x in self)

    def __sub__(self, other):
        return Weight(x - y for x, y in zip(self, other))

    def __mul__(self, other):
        return Weight(tuple(other * x for x in self))

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

    def concat(self, other):
        return Weight(self.components + other.components)
