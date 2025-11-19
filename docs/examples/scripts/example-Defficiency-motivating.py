"""
Created on Tue Jul 24 11:47:57 2012

@author: eendebakpt
"""

# %% Load necessary packages """
import sys

import matplotlib.pyplot as plt
import numpy as np

import oapackage


def addData(data, extra_data, index=0, arrays=None, extra_arrays=None):
    """Vertically stack data together"""
    nn = extra_data.shape[1]
    ww = np.vstack(
        (
            extra_data,
            index
            * np.ones(
                nn,
            ),
        )
    )
    data = np.hstack((data, ww))
    if arrays is None:
        return data
    else:
        allarrays = oapackage.joinArrayLists([arrays, extra_arrays])
        return data, allarrays


# %% Load a strength 3 orthogonal array that is D-optimal


strength3_pareto = oapackage.exampleArray(44)
arrayclass = oapackage.arraylink2arraydata(strength3_pareto)

# %% Optimize with alpha variations


def optimize_arrays(sols, alpha=None, verbose=1, niter=400):
    """Optimize set of designs"""
    vv2 = []
    ll = oapackage.arraylist_t()
    for ii, A in enumerate(sols):
        Dinit = A.Defficiency()
        Df, Ax = oapackage.optimDeffPython(A, niter=niter, verbose=0, alpha=alpha)
        D, Ds, D1 = Ax.Defficiencies()
        ll.push_back(Ax)
        oapackage.tprint("optimize array %d: %f -> D %f Ds %f D1 %f" % (ii, Dinit, D, Ds, D1))
        vv2.append(np.array([D, Ds, D1]))
    return np.array(vv2).T, ll


def getOptimizationResults(arrayclass, alpha, nrestarts=80):
    scores, dds, designs, _ = oapackage.Doptimize(arrayclass, nrestarts=nrestarts, optimfunc=alpha, selectpareto=False)
    if 1:
        dds_extra, designs_extra = optimize_arrays([strength3_pareto, strength3_pareto, strength3_pareto], alpha=alpha)
        dds = np.vstack((dds, dds_extra))
        designs = designs + list(designs_extra)
    dds = np.array([al.Defficiencies() for al in designs])
    return dds, designs


optim_functions = [[1, 0, 0], [1, 1, 0], [1, 3, 0]]

results = []
oapackage.seedfastrand(600)
oapackage.set_srand(50)
for jj, optim_func in enumerate(optim_functions):
    print("optimize series %d: optim_func %s" % (jj, optim_func))
    sys.stdout.flush()
    dds, arrays = getOptimizationResults(arrayclass, alpha=optim_func)
    results.append((dds, arrays, optim_func))


# %% Combine the data

data = np.zeros((4, 0))
allarrays = oapackage.arraylist_t()
lbls = []

for jj, result in enumerate(results):
    dds, arrays, optim_func = result
    data, allarrays = addData(data, dds.T, jj, allarrays, arrays)
    lbl = "Optimization of $D$"
    if optim_func[1] > 0:
        if optim_func[1] == 1:
            lbl += " + $D_s$"
        else:
            lbl += f" + {optim_func[1]:g}$D_s$"
    if optim_func[2] > 0:
        lbl += f" + {optim_func[2]:g}$D_{{s0}}$"
    lbls += [lbl]

data = data[[0, 1, -1], :]
x = oapackage.generateDscatter(data.T, lbls=lbls, fig=30, ndata=data.shape[0] - 1, scatterarea=80)
oapackage.oahelper.setWindowRectangle(1900, 10, w=900, h=600)

ddst3 = np.array([al.Defficiencies() for al in [strength3_pareto]])


if arrayclass.N == 40:
    plt.gca().scatter(
        ddst3[:, 1], ddst3[:, 0], marker="s", s=90, c=(0, 0, 0), linewidths=0, alpha=1, label="Strength 3"
    )
    xl = plt.xlim()
    plt.xlim([xl[0], 1.05])
    yl = plt.ylim()
    plt.ylim([0.78, yl[1]])
    plt.draw()
    plt.show()
    plt.legend(numpoints=1, scatterpoints=1, loc=3)

# researchOA.niceplot(plt.gca())
# plt.savefig(os.path.join(xdir, 'motivating-%s-scatterplot-ndata%d.png' %
#                         (arrayclass.idstr().replace('.', '-d-'), data.shape[0] - 1)), dpi=150)
