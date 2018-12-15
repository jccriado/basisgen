import itertools


# TODO: see arXiv:0909.2331
def accel_asc(n):
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        m = k + 1
        while x <= y:
            a[k] = x
            a[m] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]


def pre_partitions(n, k):
    for pre_partition in accel_asc(n):
        pre_partition = tuple(pre_partition)
        length = len(pre_partition)
        if length <= k:
            yield pre_partition + (0,) * (k - length)


def partitions(n, k):
    return set(itertools.chain.from_iterable(
        itertools.permutations(partition, k)
        for partition in pre_partitions(n, k)
    ))
