import functools


def partitions_generator(n, k):
    if k == 1:
        yield (n,)
    else:
        for m in range(n + 1):
            for previous_partition in partitions(n - m, k - 1):
                yield previous_partition + (m,)


@functools.lru_cache(maxsize=None)
def partitions(n, k):
    return list(partitions_generator(n, k))
