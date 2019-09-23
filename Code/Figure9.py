from header import *
from Time_Encoder import *
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt


def new_sig(t, delta_t, Omega):
    T = np.pi / Omega
    T_WIDTH = int(T / delta_t)
    NUM_SINCS = int(len(t) / T_WIDTH)
    # print(NUM_SINCS)
    x = np.zeros_like(t)
    sinc_loc = []
    sinc_amp = []
    for n in range(NUM_SINCS):
        sinc_loc.append(T * (n))
        a = np.random.uniform(0, 1)
        # print(a)
        sinc_amp.append(a)
        x += a * np.sinc((Omega * t - Omega * T * (n)) / np.pi) * Omega / np.pi
    scale = (np.linalg.norm(x)) / np.sqrt(len(t))
    x = x / scale
    sinc_amp = sinc_amp / scale
    return x, sinc_loc, sinc_amp


def Generate():

    kappa = 1
    delta = 1
    b = 1

    t_single = timeEncoder(kappa, delta, b, n_channels=1)
    t_half = timeEncoder(kappa, delta, b, n_channels=2, integrator_init=[-delta, 0])
    t_quarter = timeEncoder(
        kappa, delta, b, n_channels=2, integrator_init=[-delta, -0.5 * delta]
    )
    t_eighth = timeEncoder(
        kappa, delta, b, n_channels=2, integrator_init=[-delta, -0.75 * delta]
    )

    OMEGA_RANGE = [
        0.25 * np.pi,
        0.5 * np.pi,
        1 * np.pi,
        2 * np.pi,
        3 * np.pi,
        6 * np.pi,
        9 * np.pi,
        12 * np.pi,
        15 * np.pi,
    ]

    n_trials = 10

    end_time = 20
    delta_t = 1e-4
    t = np.arange(0, end_time, delta_t)
    five_percent = int(len(t) * 0.05)

    err_single = np.zeros((len(OMEGA_RANGE), n_trials))
    err_half = np.zeros((len(OMEGA_RANGE), n_trials))
    err_quarter = np.zeros((len(OMEGA_RANGE), n_trials))
    err_eighth = np.zeros((len(OMEGA_RANGE), n_trials))

    for n in range(n_trials):
        print(n)
        for o in range(len(OMEGA_RANGE)):

            Omega = OMEGA_RANGE[o]
            original_signal = new_sig(t, delta_t, Omega)[0]
            b = np.max(np.abs(original_signal)) + 1
            t_single.set_b(b)
            t_half.set_b(b)
            t_quarter.set_b(b)
            t_eighth.set_b(b)

            z_single = t_single.encode(original_signal, delta_t)
            z_half = t_half.encode(original_signal, delta_t)
            z_quarter = t_quarter.encode(original_signal, delta_t)
            z_eighth = t_eighth.encode(original_signal, delta_t)
            # print(z_half)
            # print(z_eighth)

            rec_single = t_single.decode(z_single, t, Omega, delta_t)
            rec_half = t_half.decode(z_half, t, Omega, delta_t)
            rec_quarter = t_quarter.decode(z_quarter, t, Omega, delta_t)
            rec_eighth = t_eighth.decode(z_eighth, t, Omega, delta_t)

            # print(np.linalg.norm((rec_half-rec_eighth)[five_percent:-five_percent])/(len(t)*0.9))

            err_single[o, n] = np.linalg.norm(
                (original_signal - rec_single)[five_percent:-five_percent]
            ) / (len(t) * 0.9)
            err_half[o, n] = np.linalg.norm(
                (original_signal - rec_half)[five_percent:-five_percent]
            ) / (len(t) * 0.9)
            err_quarter[o, n] = np.linalg.norm(
                (original_signal - rec_quarter)[five_percent:-five_percent]
            ) / (len(t) * 0.9)
            err_eighth[o, n] = np.linalg.norm(
                (original_signal - rec_eighth)[five_percent:-five_percent]
            ) / (len(t) * 0.9)

    # import pickle
    # with open('Different_shifts_April27.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    #    pickle.dump([error1, error2,errorN], f)

    plt.rcParams["figure.figsize"] = [8, 3]
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.weight"] = "light"
    plt.rcParams["font.size"] = "12"
    plt.rcParams["mathtext.it"] = "Times New Roman"
    plt.rcParams["mathtext.default"] = "regular"

    fig = plt.figure(figsize=(8, 3))
    plt.plot(OMEGA_RANGE, np.mean(err_single, 1), label="Single Channel")
    plt.plot(OMEGA_RANGE, np.mean(err_half, 1), label=r"shift = $\delta$")
    plt.plot(OMEGA_RANGE, np.mean(err_quarter, 1), label=r"shift = $\delta/2$")
    plt.plot(OMEGA_RANGE, np.mean(err_eighth, 1), label=r"shift = $\delta/4$")

    fig.subplots_adjust(bottom=0.2)
    ax = plt.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.xticks(OMEGA_RANGE, fontname="Times New Roman", fontsize=14)

    ax.set_xticklabels(
        [
            "0.25" r"$\pi$",
            "0.5" r"$\pi$",
            r"$\pi$",
            "2" r"$\pi$",
            "3" r"$\pi$",
            "6" r"$\pi$",
            "9" r"$\pi$",
            "12" r"$\pi$",
            "15" r"$\pi$",
        ]
    )

    plt.axvline(0.5 * np.pi, linestyle="--", color="#1f77b4")
    plt.axvline(1 * np.pi, linestyle="--", color="#ff7f0e")
    plt.text(1.2, 1e-3, "1-channel \nbound on " r"$\Omega$", rotation="vertical")

    plt.text(2.4, 1e-3, "2-channel \nbound on " r"$\Omega$", rotation="vertical")
    plt.xlabel(r"$\Omega$")
    plt.ylabel("Reconstruction error")

    plt.legend(loc="lower right")
    plt.savefig("Figures/Figure9.png")


if __name__ == "__main__":
    Generate()
