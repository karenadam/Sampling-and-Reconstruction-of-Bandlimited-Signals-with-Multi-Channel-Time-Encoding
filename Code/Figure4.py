from header import *
from Time_Encoder import *
import matplotlib
import matplotlib.pyplot as plt

from matplotlib import rc


def Generate():

    end_time = 7
    delta_t = 1e-4
    t = np.arange(0, end_time, delta_t)
    Omega = 4 * np.pi
    original_signal = (
        2 * np.cos(5 * (t + 7))
        + 0.5 * np.sin(13 * (t + 7))
        + np.cos(0.3 * (t + 7))
        + 1.2 * np.sin(0.5 * (t + 7))
        + 2 * np.cos(7 * (t + 7))
    )
    b = 6
    kappa = 1
    delta = 2
    mixing_matrix = [[1]]*2
    num_channels = 2
    t1 = timeEncoder(
        kappa,
        delta,
        b,
        mixing_matrix,
        integrator_init=[-delta, -0.25 * delta],
    )

    spikes, int_seq = t1.encode(original_signal, delta_t, with_integral_probe=True)

    if (not SimulationSettings.To_Svg):
        plt.rcParams["figure.figsize"] = [8, 6]
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["font.weight"] = "light"
        plt.rcParams["font.size"] = "12"
        plt.rcParams["mathtext.it"] = "Times New Roman"
        plt.rcParams["mathtext.default"] = "regular"

    plt.figure()
    fig, axarr = plt.subplots(4, 1)
    fig.subplots_adjust(left=0.25, right=None, hspace=0.1, wspace=-0.05)

    ax = plt.subplot(4, 1, 1)
    plt.plot(t, original_signal)
    plt.xlim(0, 7)
    ax.set_ylabel(
        "Original Signal",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    plt.gca().axes.get_xaxis().set_visible(False)

    ax = plt.subplot(4, 1, 2)
    plt.plot(t, int_seq[0, :], label=r"$y_1(t)$")
    plt.plot(t, int_seq[1, :], label=r"$y_2(t)$")
    plt.legend(loc="best")
    plt.xlim(0, 7)
    ax.set_ylabel(
        "Integrator outputs",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    plt.gca().axes.get_xaxis().set_visible(False)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax = plt.subplot(4, 1, 3)
    plt.plot(t, int_seq[0, :] - int_seq[1, :], color=colors[2])
    plt.xlim(0, 7)
    ax.set_ylabel(
        "Integrator output\n differences",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    plt.gca().axes.get_xaxis().set_visible(False)

    ax = plt.subplot(4, 1, 4)
    for ch in range(num_channels):
        for sp in spikes.get_spikes_of(ch):
            plt.axvline(sp, color=colors[ch], linestyle="--")
    plt.ylabel(
        "Spike Times",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    plt.xlabel("Time (s)")
    plt.xlim(0, 7)

    if SimulationSettings.To_Svg:
        figure_filename = Figure_Path+"Figure4_integ_comparison.svg"
        fig.savefig(figure_filename)
    else:
        figure_filename = Figure_Path+"Figure4_integ_comparison.png"
        fig.savefig(figure_filename, dpi=600)


if __name__ == "__main__":
    Generate()
