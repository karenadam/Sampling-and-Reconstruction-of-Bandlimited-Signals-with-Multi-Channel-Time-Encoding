from SimulationSettings import *
from header import *
import pickle

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

    noise_level = Figure11_noise_range[i]
    Omega = omega_range[o]
    b = np.max(np.abs(original_signal)) + 1
    tem = timeEncoder(
        kappa, delta, b, [[1]]*2, integrator_init=[-delta, 0]
    )
    sig = bandlimitedSignal(t, delta_t, Omega, sinc_locs=c1, sinc_amps=c2)
    z = tem.encode_precise(sig, Omega, end_time, tol=1e-15)
    z_corrupted = z.corrupt_with_gaussian(noise_level)
    rec = tem.decode(z_corrupted, t, Omega, delta_t, cond_n=1e-12)
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
    tem = timeEncoder(kappa, delta, b, [[1]])
    sig = bandlimitedSignal(t, delta_t, Omega, sinc_locs=c1, sinc_amps=c2)
    z = tem.encode_precise(sig, Omega, end_time, tol=1e-15)
    rec = tem.decode(z, t, Omega, delta_t, cond_n = 1e-12)
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

    err_double = np.zeros((len(omega_range), len(Figure11_noise_range), n_trials))
    err_single = np.zeros((len(omega_range), n_trials))
    five_percent = int(len(t) * 0.05)

    signals = []
    sinc_loc = []
    sinc_amp = []

    tot_number_of_runs = len(omega_range) * len(Figure11_noise_range) * n_trials
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
            get_two_channel_performance, itertools.product(signals, range(len(Figure11_noise_range)))
        )

        while not err_double.ready():
            remaining = err_double._number_left * err_double._chunksize
            sys.stderr.write("\r\033[2KRemaining: %d" % remaining)
            sys.stderr.flush()
            time.sleep(0.1)
        pool.close()
        pool.join()

    err_double = np.reshape(err_double.get(), (n_trials, len(omega_range), len(Figure11_noise_range)))
    err_double = np.transpose(err_double, (1, 2, 0))

    filename = Data_Path+"Figure11_a_VarSNR_VarBw.pkl"
    with open(filename, "wb") as f:  # Python 3: open(..., 'wb')
        pickle.dump(
            [
                err_double,
                omega_range,
                omega_range_string,
                Figure11_noise_range,
                SNR_range_string,
            ],
            f,
        )


def Generate():


    data_filename = Data_Path+"Figure11_a_VarSNR_VarBw.pkl"

    with open(data_filename, "rb") as f:  # Python 3: open(..., 'wb')
        obj = pickle.load(f, encoding="latin1")

    err_double = obj[0]
    omega_range = obj[1]
    omega_range_string = obj[2]
    Figure11_noise_range = obj[3]
    SNR_range_string = obj[4]


    data = np.mean(err_double, 2)
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
        data[::-1, :],
        ax=axes,
        norm=log_norm,
        yticklabels=omega_range_string[::-2],
        xticklabels=SNR_range_string,
        cbar_kws={"ticks": cbar_ticks[1:]},
    )
    sns_map.set_ylabel(r"Bandwidth ($\Omega$)")
    sns_map.set_xlabel("SNR")
    sns_map.set_title("Reconstruction Error")
    fig = sns_map.get_figure()
    axes = fig.gca()
    # axes.set_xticks([0, 1,2,3,4])
    plt.yticks(rotation = 0)
    axes.set_yticks([0, 2, 4, 6, 8, 10])
    fig.subplots_adjust(bottom=0.2)


    if To_Svg:
        figure_filename = Figure_Path+"Figure11_a_VarSNR_VarBw.svg"
        fig.savefig(figure_filename)
    else:
        figure_filename = Figure_Path+"Figure11_a_VarSNR_VarBw.png"
        fig.savefig(figure_filename, dpi=600)


if __name__ == "__main__":

    data_filename = Data_Path+"Figure11_a_VarSNR_VarBw.pkl"

    if not os.path.isfile(data_filename):
        GetData()
    if (SimulationSettings.graphical_import):
        Generate()
