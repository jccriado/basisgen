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
