import unittest

import numpy as np
from numpy.random import MT19937, RandomState, SeedSequence

import oapackage


def reduce(graph, colors):
    tr = oapackage.reduceGraphNauty(graph, colors=colors, verbose=0)
    tri = inverse_permutation(tr)
    graph_reduced = oapackage.transformGraphMatrix(graph, tri)
    colors_reduced = [colors[idx] for idx in tr]
    return graph_reduced, colors_reduced, tr


def inverse_permutation(perm):
    inverse = [0] * len(perm)
    for i, p in enumerate(perm):
        inverse[p] = i
    return inverse


def test_reduceNauty():
    rs = RandomState(MT19937(SeedSequence(123456789)))

    for graph_size in [2, 4, 8]:
        graph = (rs.random((graph_size, graph_size)) > 0.5).astype(int)
        colors = [int(v) for v in rs.random(graph_size) > .5]

        graph_reduced, colors_reduced, tr = reduce(graph, colors)
        perm = rs.permutation(graph_size)

        graph2 = graph[perm, :][:, perm]
        colors2 = [colors[idx] for idx in perm]

        graph2_reduced, colors2_reduced, tr2 = reduce(graph2, colors2)

        if not np.all(graph_reduced == graph2_reduced):
            raise Exception('reduced graphs unequal')
        if not np.all(colors_reduced == colors2_reduced):
            raise Exception('reduced colors unequal')


if __name__ == '__main__':
    """ Test code """
    unittest.main()
