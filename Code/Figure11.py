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
    s = args[0][1]
    original_signal = args[0][2]
    z = args[0][3]
    kappa = args[0][4]
    delta = args[0][5]
    end_time = args[0][6]
    delta_t = args[0][7]
    t = np.arange(0, end_time, delta_t)
    five_percent = int(len(t) * 0.05)
    l = args[1]

    shift = shifts[s]
    noise_level = noise_range[l]

    z_corrupted = z.corrupt_with_gaussian(noise_level)
    b = np.max(np.abs(original_signal)) + 1
    tem = timeEncoder(kappa, delta, b, [[1]]*2,  integrator_init=[-delta, -delta + shift * delta])

    # print(z_corrupted.get_spikes_of(0))
    # print(z.get_spikes_of(0))
    rec = tem.decode(z_corrupted, t, Omega_var_noise, delta_t)
    err = np.linalg.norm((original_signal - rec)[five_percent:-five_percent]) / (
        len(t) * 0.9
    )
    return err


def get_single_channel_performance(args):

    n = args[0][0]
    original_signal = args[0][1]
    z_single = args[0][2]
    kappa = args[0][3]
    delta = args[0][4]
    end_time = args[0][5]
    delta_t = args[0][6]
    t = np.arange(0, end_time, delta_t)
    five_percent = int(len(t) * 0.05)
    l = args[1]

    noise_level = noise_range[l]

    z_corrupted = z_single.corrupt_with_gaussian(noise_level)
    b = np.max(np.abs(original_signal)) + 1
    tem = timeEncoder(kappa, delta, b, [[1]])

    # print(z_corrupted.get_spikes_of(0))
    # print(z.get_spikes_of(0))
    rec = tem.decode(z_corrupted, t, Omega_var_noise, delta_t)
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

    signals_multi_channel = []
    signals_single_channel = []

    tot_number_of_runs = len(omega_range) * len(shifts) * n_trials
    run_tracker_divisor = int(tot_number_of_runs / 25) + 1

    for n in range(n_trials):
        original_signal, c1, c2 = new_sig(t, delta_t, Omega_var_noise)
        b = np.max(np.abs(original_signal)) + 1

        for s in range(len(shifts)):
            shift = shifts[s]
            tem = timeEncoder(kappa, delta, b, [[1]]*2, integrator_init=[-delta, -delta + shift * delta])
            sig = bandlimitedSignal(t, delta_t, Omega_var_noise, sinc_locs=c1, sinc_amps=c2)
            z = tem.encode_precise(sig, Omega_var_noise, end_time, tol=1e-15)
            signals_multi_channel.append([n,s,original_signal, z, kappa, delta, end_time, delta_t])
           
        z_single = z.get_spikes_of(0, asSpikeTimesObject = True)  
        signals_single_channel.append([n,original_signal,z_single, kappa, delta, end_time, delta_t])

    with multiprocessing.Pool(processes=int(multiprocessing.cpu_count() / 2)) as pool:
        err_double = pool.map_async(
            get_two_channel_performance, itertools.product(signals_multi_channel, range(len(noise_range)))
        )
        err_single = pool.map_async(
            get_single_channel_performance, itertools.product(signals_single_channel, range(len(noise_range)))
        )
        while not err_double.ready():
            remaining = err_double._number_left * err_double._chunksize
            sys.stderr.write("\r\033[2KRemaining: %d" % remaining)
            sys.stderr.flush()
            time.sleep(0.1)
        pool.close()
        pool.join()

    err_double = np.reshape(err_double.get(), (n_trials, len(shifts), len(noise_range)))
    err_double = np.transpose(err_double, (2, 1, 0))
    err_single = np.reshape(err_single.get(), (n_trials, len(noise_range)))
    err_single = np.transpose(err_single)

    filename = Data_Path+"Figure11_VarShifts.pkl"
    with open(filename, "wb") as f:  # Python 3: open(..., 'wb')
        pickle.dump(
            [
                err_double,
                err_single,
                noise_range,
                noise_range_string,
                shifts,
                shifts_range_string,
            ],
            f,
        )


def Generate():


    data_filename = Data_Path+"Figure11_VarShifts.pkl"

    with open(data_filename, "rb") as f:  # Python 3: open(..., 'wb')
        obj = pickle.load(f, encoding="latin1")

    err_double = obj[0]
    err_single = obj[1]
    noise_range = obj[2]
    noise_range_string = obj[3]
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
        yticklabels=SNR_range_string[::-1],
        xticklabels=shifts_range_string[::2],
        cbar_kws={"ticks": cbar_ticks[1:]},
    )
    sns_map.set_ylabel(r"SNR")
    sns_map.set_xlabel("Shift")
    sns_map.set_title("Reconstruction Error")
    fig = sns_map.get_figure()
    axes = fig.gca()
    axes.set_xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
    axes.set_ylim(6,0)
    axes.set_yticklabels(SNR_range_string[::-1], rotation = 0)
    # axes.set_yticks([0, 2, 4, 6, 8, 10, 12])
    axes.vlines(17, *axes.get_ylim(), color="yellow", linestyle="--")
    fig.subplots_adjust(bottom=0.2)
    if To_Svg:
        figure_filename = Figure_Path+"Figure11_VarShifts.svg"
        fig.savefig(figure_filename)
    else:
        figure_filename = Figure_Path+"Figure11_VarShifts.png"
        fig.savefig(figure_filename, dpi=600)


if __name__ == "__main__":

    data_filename = Data_Path+"Figure11_VarShifts.pkl"

    if not os.path.isfile(data_filename):
        GetData()
    if (SimulationSettings.graphical_import):
        Generate()
