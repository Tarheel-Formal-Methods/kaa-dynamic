from kaa.experiutil import generate_init_traj
from sklearn.decomposition import PCA
from models.sir import SIR

import numpy as np
import matplotlib.pyplot as plt

def test_pca():

    model = SIR()
    trajs = generate_init_traj(model, 20, 40)

    traj_mat = np.empty((20,3))

    for traj_idx, traj in enumerate(trajs):
        traj_mat[traj_idx] = traj[-1]

    pca = PCA(n_components=3)
    pca.fit(traj_mat)

    print("Components: {} \n Mean: {}".format(pca.components_, pca.mean_))

    # plot data
    #plt.scatter(traj_mat[:, 0], traj_mat[:, 1], alpha=0.2)
    #for length, vector in zip(pca.explained_variance_, pca.components_):
    #    v = vector * 3 * np.sqrt(length)
    #    draw_vector(pca.mean_, pca.mean_ + v)
    #    plt.axis('equal')
    #plt.show()


def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops=dict(arrowstyle='->',
                    linewidth=2,
                    shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)
