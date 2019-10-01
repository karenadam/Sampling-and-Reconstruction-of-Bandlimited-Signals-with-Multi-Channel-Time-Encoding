from SimulationSettings import *
from header import *
import pickle
from matplotlib.colors import LogNorm

try:
    from ExtraSpecifications import *
except:
    pass


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


def get_two_channel_performance(args):

    n = args[0][0]
    o = args[0][1]
    original_signal = args[0][2]
    c1 = args[0][3]
    c2 = args[0][4]
    kappa = args[0][5]
    delta = args[0][6]
    end_time = args[0][7]
    delta_t = args[0][8]
    t = np.arange(0, end_time, delta_t)
    five_percent = int(len(t) * 0.05)
    i = args[1]

    shift = shifts[i]
    Omega = omega_range[o]
    b = np.max(np.abs(original_signal)) + 1
    tem = timeEncoder(
        kappa, delta, b, n_channels=2, integrator_init=[-delta, -delta + shift * delta]
    )
    z = tem.encode_precise(c1, c2, Omega, end_time)
    rec = tem.decode(z, t, Omega, delta_t, tol=1e-15)
    err = np.linalg.norm((original_signal - rec)[five_percent:-five_percent]) / (
        len(t) * 0.9
    )
    return err


def get_single_channel_performance(args):

    n = args[0][0]
    o = args[0][1]
    original_signal = args[0][2]
    c1 = args[0][3]
    c2 = args[0][4]
    kappa = args[0][5]
    delta = args[0][6]
    end_time = args[0][7]
    delta_t = args[0][8]
    t = np.arange(0, end_time, delta_t)
    five_percent = int(len(t) * 0.05)

    Omega = omega_range[o]
    b = np.max(np.abs(original_signal)) + 1
    tem = timeEncoder(kappa, delta, b, n_channels=1)
    z = tem.encode_precise(c1, c2, Omega, end_time)
    rec = tem.decode(z, t, Omega, delta_t, tol=1e-15)
    err = np.linalg.norm((original_signal - rec)[five_percent:-five_percent]) / (
        len(t) * 0.9
    )
    return err


def GetData():

    np.random.seed(0)

    kappa = 1
    delta = 1

    end_time = 20
    delta_t = 1e-4
    t = np.arange(0, end_time, delta_t)

    err_double = np.zeros((len(omega_range), len(shifts), n_trials))
    err_single = np.zeros((len(omega_range), n_trials))
    five_percent = int(len(t) * 0.05)

    signals = []
    sinc_loc = []
    sinc_amp = []

    tot_number_of_runs = len(omega_range) * len(shifts) * n_trials
    run_tracker_divisor = int(tot_number_of_runs / 25) + 1

    for n in range(n_trials):
        for o in range(len(omega_range)):
            Omega = omega_range[o]
            original_signal, c1, c2 = new_sig(t, delta_t, Omega)
            signals.append(
                [n, o, original_signal, c1, c2, kappa, delta, end_time, delta_t]
            )
            sinc_loc.append(c1)
            sinc_amp.append(c2)

    with multiprocessing.Pool(processes=int(multiprocessing.cpu_count() / 2)) as pool:
        err_double = pool.map_async(
            get_two_channel_performance, itertools.product(signals, range(len(shifts)))
        )
        err_single = pool.map_async(
            get_single_channel_performance, itertools.product(signals, repeat=1)
        )
        while not err_double.ready():
            remaining = err_double._number_left * err_double._chunksize
            sys.stderr.write("\r\033[2KRemaining: %d" % remaining)
            sys.stderr.flush()
            time.sleep(0.1)
        pool.close()
        pool.join()

    err_double = np.reshape(err_double.get(), (n_trials, len(omega_range), len(shifts)))
    err_double = np.transpose(err_double, (1, 2, 0))
    err_single = np.reshape(err_single.get(), (n_trials, len(omega_range)))
    err_single = np.transpose(err_single)

    filename = Data_Path+"Figure10_VarShifts.pkl"
    with open(filename, "wb") as f:  # Python 3: open(..., 'wb')
        pickle.dump(
            [
                err_double,
                err_single,
                omega_range,
                omega_range_string,
                shifts,
                shifts_range_string,
            ],
            f,
        )


def Generate():

    figure_filename = Figure_Path+"Figure10_VarShifts.png"

    data_filename = Data_Path+"Figure10_VarShifts.pkl"

    with open(data_filename, "rb") as f:  # Python 3: open(..., 'wb')
        obj = pickle.load(f, encoding="latin1")

    err_double = obj[0]
    err_single = obj[1]
    omega_range = obj[2]
    omega_range_string = obj[3]
    shifts = obj[4]
    shifts_range_string = obj[5]

    shifts_range_string.append(0)
    shifts_range_string.append(0)

    data = np.mean(err_double, 2)
    data_aug = np.zeros((data.shape[0], data.shape[1] + 2))
    data_aug[:, : data.shape[1]] = data[:, :]
    data_aug[:, -2] = np.mean(err_single, 1)
    data_aug[:, -1] = np.mean(err_single, 1)
    log_norm = LogNorm(vmin=data.min().min(), vmax=data.max().max())
    cbar_ticks = [
        math.pow(10, i)
        for i in range(
            math.floor(math.log10(data.min().min())),
            1 + math.ceil(math.log10(data.max().max())),
        )
    ]

    fig = plt.figure(figsize=(4, 3.3))
    fig.subplots_adjust(left=0.15)
    axes = fig.gca()
    sns_map = sns.heatmap(
        data_aug[::-1, :],
        ax=axes,
        norm=log_norm,
        yticklabels=omega_range_string[::-2],
        xticklabels=shifts_range_string[::2],
        cbar_kws={"ticks": cbar_ticks[1:]},
    )
    sns_map.set_ylabel(r"Bandwidth ($\Omega$)")
    sns_map.set_xlabel("Shift")
    sns_map.set_title("Reconstruction Error")
    fig = sns_map.get_figure()
    axes = fig.gca()
    axes.set_xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22])
    axes.set_yticks([0, 2, 4, 6, 8, 10, 12])
    axes.vlines(21, *axes.get_ylim(), color="yellow", linestyle="--")
    fig.subplots_adjust(bottom=0.2)
    fig.savefig(figure_filename, dpi=600)


if __name__ == "__main__":

    data_filename = Data_Path+"Figure10_VarShifts.pkl"

    if not os.path.isfile(data_filename):
        GetData()
    Generate()
