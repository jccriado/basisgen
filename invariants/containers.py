import collections


class MultivaluedMap(collections.defaultdict):
    def __init__(self, *args, **kwargs):
        return super().__init__(set, *args, **kwargs)

    @staticmethod
    def from_pairs(pairs):
        mapping = MultivaluedMap()
        for key, value in pairs:
            mapping[key].add(value)

        return mapping

    def update(self, other):
        for key, value in other.items():
            self[key].update(value)


class OrderedCounter(collections.Counter, collections.OrderedDict):
    pass


class Tree(object):
    def __init__(self, node, left, right):
        self.node = node
        self.left = left
        self.right = right

    def __str__(self):
        left = [line for line in str(self.left).split('\n') if line != '']
        head_left = "|-- {}".format(left[0])
        tail_left = "\n".join("|   {}".format(line) for line in left[1:])

        right = [line for line in str(self.right).split('\n') if line != '']
        head_right = "`-- {}".format(right[0])
        tail_right = "\n".join("    {}".format(line) for line in right[1:])

        return (
            f"{self.node}"
            + f"\n{head_left}\n{tail_left}\n{head_right}\n{tail_right}"
        )

    __repr__ = __str__
