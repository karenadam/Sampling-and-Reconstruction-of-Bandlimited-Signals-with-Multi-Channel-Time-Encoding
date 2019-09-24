from header import *
from Time_Encoder import *
import matplotlib
import matplotlib.pyplot as plt

from brian2 import *
prefs.codegen.target = 'numpy' 
from matplotlib import rc


def Generate():
    matplotlib.rc("text", usetex=True)
    matplotlib.rc("font", family="serif")
    matplotlib.rc("font", size=24)
    matplotlib.rc("text.latex", preamble=r"\usepackage{amsmath}\usepackage{amssymb}")

    #####  ORIGINAL INPUT SIGNAL  #####
    ###################################
    delta_t = 1e-4
    time = np.arange(0, 10, delta_t)
    signal = 2.5 + np.cos(0.3 * np.pi * time - 2) - np.sin(0.1 * np.pi * time - 1)

    sampled_time = time[1 : len(signal) : 5000]
    sampled_signal = signal[1 : len(signal) : 5000]

    plt.figure()
    plt.plot(time, signal)
    plt.ylim(0, 4)
    plt.xlim(0, 10)
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.savefig(Figure_Path+"Figure1_orig_signal.svg")

    #####   CLASSICALLY SAMPLED   #####
    ###################################

    plt.figure()
    plt.stem(sampled_time, sampled_signal, use_line_collection = True)
    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.ylim(0, 4)
    plt.xlim(0, 10)
    plt.savefig(Figure_Path+"Figure1_class_samples.svg")

    #####       TIME ENCODED      #####
    ###################################

    t = timeEncoder(1, 0.75, 1, 1, [0])
    z = t.encode(signal, delta_t)
    spikes = z.get_spikes_of(0)
    plt.figure()
    plt.stem(spikes, np.ones_like(spikes), use_line_collection = True)
    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xlabel("Time (s)")
    plt.ylabel(" ")
    plt.savefig(Figure_Path+"Figure1_time_encoding.svg")

    #####    PASSED THROUGH LIF   #####
    ###################################
    @check_units(v=volt, tol=1, result=volt / (second ** 2))
    def dirac(v):
        if (v >= -1e-4 * volt) and (v <= 1e-4 * volt):
            return 1 * volt / (second ** 2)
        else:
            return 0 * volt / (second ** 2)

    n = 1
    duration = 1 * second
    tau = 100 * ms
    Cm = 20 * uF
    Vreset = 0 * volt
    Vth = 10e-3 * volt
    eqs = """
    Iinj = 2.5e-6*amp + A1*cos(freq1*pi*t-2) - A2*sin(freq2*pi*t-1) : amp (constant over dt)
    Ileak = -Cm/tau *(v-v0) : amp 
    Ispike = Cm* (tau/(v0-Vth))*(Vreset-Vth)*dirac(v-Vth): amp (constant over dt)
    dv/dt = 1/Cm*(Ileak+Iinj) : volt (unless refractory)
    v0 : volt
    freq1 : hertz/radian
    freq2 : hertz/radian
    A1 : amp
    A2 : amp
    """
    group = NeuronGroup(
        n,
        eqs,
        threshold="v > 10*mV",
        reset="v = 0*mV",
        refractory=5 * ms,
        method="exact",
    )
    group.v = 0 * mV
    group.v0 = "20*mV"

    group.freq1 = 3 / second
    group.freq2 = 1 / second
    group.A1 = 1e-6 * amp
    group.A2 = 1e-6 * amp
    group.Iinj = 1.5 * (1.5e-5 * amp + group.A1 * cos(-2) - group.A2 * sin(-1))

    monitor = SpikeMonitor(group)
    mon = StateMonitor(group, ["v", "Ispike", "Iinj", "Ileak"], record=True)

    run(duration)
    plt.figure()
    plt.plot(mon.t, (-mon.Ileak[0] - mon.Ispike[0]) * 1e6)
    plt.ylabel(r"Neuron Output Current ($\mu A$)")
    plt.xlabel("Time (s)")
    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.savefig(Figure_Path+"Figure1_LIF_Sampling.svg")



if __name__ == "__main__":
    Generate()
