"""Example script demonstrating D-efficiency optimization with different alpha values.

This script shows how different optimization objectives (D, Ds, D1) affect the
resulting designs when optimizing orthogonal arrays.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

import oapackage


def add_data(
    data: np.ndarray,
    extra_data: np.ndarray,
    index: int = 0,
    arrays: oapackage.arraylist_t | None = None,
    extra_arrays: list | None = None,
) -> np.ndarray | tuple[np.ndarray, oapackage.arraylist_t]:
    """Vertically stack data together with an index row.

    Args:
        data: Existing data array to append to.
        extra_data: New data to add.
        index: Index value to add as a row.
        arrays: Optional existing array list.
        extra_arrays: Optional new arrays to join.

    Returns:
        Updated data array, or tuple of (data, combined arrays) if arrays provided.
    """
    nn = extra_data.shape[1]
    ww = np.vstack((extra_data, index * np.ones(nn)))
    data = np.hstack((data, ww))

    if arrays is None:
        return data

    all_arrays = oapackage.joinArrayLists([arrays, extra_arrays])
    return data, all_arrays


# %% Load a strength 3 orthogonal array that is D-optimal


strength3_pareto = oapackage.exampleArray(44)
arrayclass = oapackage.arraylink2arraydata(strength3_pareto)

# %% Optimize with alpha variations


def optimize_arrays(
    sols: list, alpha: list | None = None, niter: int = 400
) -> tuple[np.ndarray, oapackage.arraylist_t]:
    """Optimize a set of designs.

    Args:
        sols: List of arrays to optimize.
        alpha: Alpha weights for optimization objective.
        niter: Number of optimization iterations.

    Returns:
        Tuple of (efficiencies array, optimized array list).
    """
    efficiencies = []
    optimized_arrays = oapackage.arraylist_t()

    for ii, array in enumerate(sols):
        d_init = array.Defficiency()
        _, optimized = oapackage.optimDeffPython(array, niter=niter, verbose=0, alpha=alpha)
        d, ds, d1 = optimized.Defficiencies()
        optimized_arrays.push_back(optimized)
        oapackage.tprint(f"optimize array {ii}: {d_init:f} -> D {d:f} Ds {ds:f} D1 {d1:f}")
        efficiencies.append(np.array([d, ds, d1]))

    return np.array(efficiencies).T, optimized_arrays


def get_optimization_results(
    arrayclass: oapackage.arraydata_t, alpha: list, nrestarts: int = 80
) -> tuple[np.ndarray, list]:
    """Run D-optimization and return results.

    Args:
        arrayclass: Array class specification.
        alpha: Alpha weights for optimization objective.
        nrestarts: Number of random restarts.

    Returns:
        Tuple of (D-efficiencies array, list of designs).
    """
    r = oapackage.Doptimize(arrayclass, nrestarts=nrestarts, alpha=alpha, verbose=False)
    designs = r.designs
    _, designs_extra = optimize_arrays([strength3_pareto, strength3_pareto, strength3_pareto], alpha=alpha)
    designs = list(designs) + list(designs_extra)
    dds = np.array([al.Defficiencies() for al in designs])

    return dds, designs


optim_functions = [[1, 0, 0], [1, 1, 0], [1, 3, 0]]

oapackage.seedfastrand(600)
oapackage.set_srand(50)
results = []
for jj, optim_func in enumerate(optim_functions):
    print(f"optimize series {jj}: optim_func {optim_func}")
    dds, arrays = get_optimization_results(arrayclass, alpha=optim_func)
    results.append((dds, arrays, optim_func))


# %% Combine the data

data = np.zeros((4, 0))
all_arrays = oapackage.arraylist_t()
labels = []

for jj, result in enumerate(results):
    dds, arrays, optim_func = result
    data, all_arrays = add_data(data, dds.T, jj, all_arrays, arrays)
    label = "Optimization of $D$"
    if optim_func[1] > 0:
        if optim_func[1] == 1:
            label += " + $D_s$"
        else:
            label += f" + {optim_func[1]:g}$D_s$"
    if optim_func[2] > 0:
        label += f" + {optim_func[2]:g}$D_{{s0}}$"
    labels.append(label)

data = data[[0, 1, -1], :]
oapackage.generateDscatter(data.T, lbls=labels, fig=30, ndata=data.shape[0] - 1, scatterarea=80)
oapackage.oahelper.setWindowRectangle(1900, 10, w=900, h=600)

strength3_efficiencies = np.array([strength3_pareto.Defficiencies()])


if arrayclass.N == 40:
    plt.gca().scatter(
        strength3_efficiencies[:, 1],
        strength3_efficiencies[:, 0],
        marker="s",
        s=90,
        color=(0, 0, 0),
        linewidths=0,
        alpha=1,
        label="Strength 3",
    )
    xl = plt.xlim()
    plt.xlim([xl[0], 1.05])
    yl = plt.ylim()
    plt.ylim([0.78, yl[1]])
    plt.draw()
    plt.show()
    plt.legend(numpoints=1, scatterpoints=1, loc=3)

oapackage.niceplot(plt.gca())
