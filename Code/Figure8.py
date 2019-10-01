from header import *
from SimulationSettings import *
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


def get_M_channel_performance(args):
    kappa = 1
    delta = 1

    err_mult = np.zeros((len(omega_range_M_channels), len(num_channels), n_trials))

    n = args[0][0]
    o = args[0][1]
    original_signal = args[0][2]
    c1 = args[0][3]
    c2 = args[0][4]
    i = args[1]

    five_percent = int(len(original_signal) * 0.05)
    t = np.arange(0, end_time, delta_t)

    M = num_channels[i]
    Omega = omega_range_M_channels[o]
    b = np.max(np.abs(original_signal)) + 1
    integ_init = np.arange(-delta, delta, 2 * delta / M).tolist()

    tem2 = timeEncoder(kappa, delta, b, n_channels=M, integrator_init=integ_init)
    spikes = tem2.encode(original_signal, delta_t)
    rec = tem2.decode(spikes, t, Omega, delta_t)
    # rec = tem.decode(z2, t, Omega, delta_t)
    err = np.linalg.norm((original_signal - rec)[five_percent:-five_percent]) / (
        len(t) * 0.9
    )
    return err


def GetData():
    np.random.seed(0)
    signals = []
    sinc_loc = []
    sinc_amp = []

    tot_number_of_runs = len(omega_range_M_channels) * len(shifts) * n_trials
    run_tracker_divisor = int(tot_number_of_runs / 25) + 1

    end_time = 20
    delta_t = 1e-4
    t = np.arange(0, end_time, delta_t)

    for n in range(n_trials):
        for o in range(len(omega_range_M_channels)):
            Omega = omega_range_M_channels[o]
            original_signal, c1, c2 = new_sig(t, delta_t, Omega)
            signals.append([n, o, original_signal, c1, c2])
            sinc_loc.append(c1)
            sinc_amp.append(c2)

    with multiprocessing.Pool(processes=int(multiprocessing.cpu_count() / 2)) as pool:
        err_mult = pool.map_async(
            get_M_channel_performance,
            itertools.product(signals, range(len(num_channels))),
        )
        while not err_mult.ready():
            remaining = err_mult._number_left * err_mult._chunksize
            sys.stderr.write("\r\033[2KRemaining: %d" % remaining)
            sys.stderr.flush()
            time.sleep(0.1)
        pool.close()
        pool.join()

    err_mult = np.reshape(
        err_mult.get(), (n_trials, len(omega_range_M_channels), len(num_channels))
    )
    err_mult = np.transpose(err_mult, (1, 2, 0))

    import pickle

    filename = Data_Path+"Figure8_VarNumChannels.pkl"
    
    with open(filename, "wb") as f:  # Python 3: open(..., 'wb')
        pickle.dump(
            [
                err_mult,
                omega_range_M_channels,
                omega_range_M_channels_string,
                num_channels,
            ],
            f,
        )


def Generate():
    figure_filename = Figure_Path+"Figure8_VarNumChannels.png"

    data_filename = Data_Path+"Figure8_VarNumChannels.pkl"

    with open(data_filename, 'rb') as f:  # Python 3: open(..., 'wb')
        obj = pickle.load(f,encoding='latin1')

    err_Mult = obj[0]
    omega_range = obj[1]
    omega_range_string = obj[2] 
    number_of_channels = obj[3]

 
        
    data = np.mean(err_Mult,2)
    log_norm = LogNorm(vmin=data.min().min(), vmax=data.max().max())
    cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(data.min().min())), 1+math.ceil(math.log10(data.max().max())))]

    fig = plt.figure(figsize = (4,3.3))
    fig.subplots_adjust(left = 0.15)
    axes = fig.gca()
    sns_map = sns.heatmap(data[::-1,:],ax = axes, norm=log_norm,yticklabels = omega_range_string[::-2],xticklabels = number_of_channels, cbar_kws = {'ticks':cbar_ticks})
    sns_map.set_ylabel(r'Bandwidth ($\Omega$)')
    sns_map.set_xlabel('Number of Channels')
    sns_map.set_title('Reconstruction Error')
    plt.yticks(rotation = 0)
    axes.set_yticks(np.arange(0,len(omega_range),2))
    axes.vlines(21, *axes.get_ylim(), color = 'yellow', linestyle = '--')
    fig.subplots_adjust(bottom = 0.2)

    if To_Svg:
        figure_filename = Figure_Path+"Figure8_VarNumChannels.svg"
        fig.savefig(figure_filename)
    else:
        figure_filename = Figure_Path+"Figure8_VarNumChannels.png"
        fig.savefig(figure_filename, dpi=600)




if __name__ == "__main__":

    data_filename = "Data/Figure8_VarNumChannels.pkl"

    if not os.path.isfile(data_filename):
        GetData()
    if (SimulationSettings.graphical_import):
        Generate()
