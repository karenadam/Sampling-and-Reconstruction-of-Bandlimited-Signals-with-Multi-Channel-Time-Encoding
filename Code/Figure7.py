from header import *
from Time_Encoder import *
import matplotlib
import matplotlib.pyplot as plt

from brian2 import *
from matplotlib import rc


def new_signal(Omega, num_freq=0):
    """ Creates a new signal with bandwidth Omega"""
    """ Signal has at most 20 cosines, frequency content is assigned randomly"""
    """ Signal energy is normalized to 1"""
    if num_freq > 0:
        NUM_COS = num_freq
    else:
        NUM_COS = np.random.randint(0, 50) + 1
    cosfreq = np.zeros((NUM_COS, 1))
    cosamp = np.zeros((NUM_COS, 1))
    for i in range(NUM_COS):
        cosfreq[i] = np.random.uniform(0, Omega)
        cosamp[i] = np.random.uniform(0, 1)

    cosfreq = cosfreq * Omega / max(cosfreq)
    cosamp = cosamp / (np.linalg.norm(np.transpose(cosamp)))

    return cosfreq, cosamp


def sample_signal(t, cosfreq, cosamp, b=0):
    signal = np.zeros_like(t)
    for i in range(len(cosfreq)):
        signal += cosamp[i] * (np.cos(cosfreq[i] * t)) + b
    return signal


def Generate():
    kappa = 4
    delta = 1

    omega = np.pi
    delta_t = 1e-2
    t = np.arange(0, 25, delta_t)
    np.random.seed(10)
    c1, c2 = new_signal(omega)
    y = sample_signal(t, c1, c2) + 1
    b = np.max(np.abs(y)) + 1

    tem_double = timeEncoder(kappa, delta, b, [[1]]*2, integrator_init=[-delta, 0])
    spikes_double = tem_double.encode(y, delta_t)

    rec_dbl = tem_double.decode(spikes_double, t, omega, delta_t)

    spikes_single = spikes_double.get_spikes_of(0, asSpikeTimesObject=True)
    tem_single = timeEncoder(kappa, delta, b, [1])
    rec2 = tem_single.decode(spikes_single, t, omega, delta_t)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    plt.figure()
    fig, axarr = plt.subplots(4, 1)

    fig.subplots_adjust(left=0.25, right=None, hspace=0.1, wspace=-0.05)

    ax = plt.subplot(3, 1, 1)
    plt.plot(t, y, linewidth=1, label="Original Signal", color=colors[0])
    plt.plot(
        t, rec2, linewidth=1, linestyle="--", label="Reconstruction", color=colors[2]
    )
    plt.xlim(5, 20)
    plt.ylim(-2, 5)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.legend(loc="upper left")
    ax.set_ylabel(
        "Reconstruction\n using one TEM",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    ax.yaxis.set_label_coords(-0.1, 0.5)

    ax = plt.subplot(3, 1, 2)
    plt.plot(t, y, linewidth=1, label="Original Signal", color=colors[0])
    plt.plot(
        t, rec_dbl, linewidth=1, linestyle="--", label="Reconstruction", color=colors[1]
    )
    plt.xlim(5, 20)
    plt.ylim(-2, 5)
    plt.legend(loc="upper left")
    ax.set_ylabel(
        "Reconstruction\n   using two TEMs",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    plt.xlabel("Time (s)")

    ax = plt.subplot(3, 1, 3)
    plt.plot(
        t,
        rec2 - y,
        linewidth=1,
        linestyle="--",
        label="1-Channel Reconstruction Error",
        color=colors[2],
    )
    plt.plot(
        t,
        rec_dbl - y,
        linewidth=1,
        linestyle="--",
        label="2-Channel Reconstruction Error",
        color=colors[1],
    )
    plt.xlim(5, 20)
    plt.ylim(-2, 5)
    plt.legend(loc="upper left")
    ax.set_ylabel(
        "Reconstruction\n   error",
        rotation=-1,
        horizontalalignment="right",
        verticalalignment="center",
    )
    plt.xlabel("Time (s)")
    
    if SimulationSettings.To_Svg:
        figure_filename = Figure_Path+"Figure7_MChannelRec.svg"
        fig.savefig(figure_filename)
    else:
        figure_filename = Figure_Path+"Figure7_MChannelRec.png"
        fig.savefig(figure_filename, dpi=600)


if __name__ == "__main__":
    Generate()
