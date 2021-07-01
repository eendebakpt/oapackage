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


class TestGWLP(unittest.TestCase):

    def test_reduceNauty_directed_graphs(self):
        rs = RandomState(MT19937(SeedSequence(123456789)))

        for jj, graph_size in enumerate([2, 4, 5, 6, 8]):
            graph = (rs.random((graph_size, graph_size)) > 0.5).astype(int)
            colors = [int(v) for v in rs.random(graph_size) > .5]

            graph_reduced, colors_reduced, tr = reduce(graph, colors)

            for ii in range(4):
                perm = rs.permutation(graph_size)

                graph2 = graph[perm, :][:, perm]
                colors2 = [colors[idx] for idx in perm]
                graph2_reduced, colors2_reduced, tr2 = reduce(graph2, colors2)

                self.assertTrue(np.all(graph_reduced == graph2_reduced))
                self.assertTrue(np.all(colors_reduced == colors2_reduced))


if __name__ == '__main__':
    """ Test code """
    unittest.main()
