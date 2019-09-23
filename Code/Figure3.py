from header import *
from Time_Encoder import *
import matplotlib
import matplotlib.pyplot as plt

from brian2 import *
from matplotlib import rc


def Generate():

    matplotlib.rc("text", usetex=True)
    matplotlib.rc("font", family="serif")
    matplotlib.rc("font", size=12)
    matplotlib.rc("text.latex", preamble=r"\usepackage{amsmath}\usepackage{amssymb}")

    # matplotlib.rcParams['mathtext.fontset'] = 'stix'
    delta_t = 1e-2
    t = np.arange(0, 10, delta_t)
    x = (
        2 * np.cos(5 * t)
        + 0.5 * np.sin(13 * t)
        + np.cos(0.3 * t)
        + 1.2 * np.sin(0.5 * t)
        + 2 * np.cos(7 * t)
    )

    b = 6
    kappa = 1
    delta = 2

    x_b = x + b
    int_x_b = (np.cumsum(x_b) * delta_t / kappa) % (2 * delta) - delta
    tem = timeEncoder(kappa, delta, b)
    z = tem.encode(x, delta_t)
    spikes = z.get_spikes_of(0)
    plt.figure()
    plt.tight_layout()

    fig, axarr = plt.subplots(4, 1)

    fig.subplots_adjust(left=0.25, right=None, hspace=0.1, wspace=-0.05)

    ax = plt.subplot(4, 1, 1)
    plt.plot(t, x, linewidth=1)
    plt.xlim(0, 10)
    plt.gca().axes.get_xaxis().set_visible(False)
    ax.set_ylabel(
        r"$x(t)$", rotation=-1, horizontalalignment="right", verticalalignment="center"
    )
    plt.ylim(-7, 13)
    ax.yaxis.set_label_coords(-0.1, 0.5)

    ax = plt.subplot(4, 1, 2)
    plt.plot(t, x_b, linewidth=1)
    plt.xlim(0, 10)
    plt.gca().axes.get_xaxis().set_visible(False)
    ax.set_ylabel(
        r"$x(t)+b$",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    plt.ylim(-7, 13)
    ax.yaxis.set_label_coords(-0.1, 0.5)

    ax = plt.subplot(4, 1, 3)
    p = plt.plot(t, int_x_b, linewidth=1)
    color = ax.get_lines()[0].get_c()
    plt.xlim(0, 10)
    plt.gca().axes.get_xaxis().set_visible(False)
    ax.set_ylabel(
        r"$\int_{t_k}^t (x(u)+b)\, du$",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    ax.yaxis.set_label_coords(-0.1, 0.5)

    ax = plt.subplot(4, 1, 4)
    for spike in spikes:
        plt.axvline(spike, linewidth=1, linestyle="--")
    # plt.stem(z,10*np.ones_like(z), linefmt=color+'--')
    plt.xlim(0, 10)
    plt.ylim(0.1, 1)
    # plt.gca().axes.get_yaxis().set_visible(False)
    ax.set_ylabel(
        r"$\left\{t_k, k \in \mathbb{Z}\right\}$",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    ax.yaxis.set_label_coords(-0.1, 0.5)
    plt.yticks([])
    plt.xlabel("Time (s)")

    plt.savefig("../Figures/Figure3_TEExample.jpg", dpi=300)


if __name__ == "__main__":
    Generate()
