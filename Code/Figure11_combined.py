from SimulationSettings import *
# from  SimulationSettings import Figure11_Figure11_Omega_Range as Figure11_Omega_Range
# from  SimulationSettings import Figure11_Figure11_Omega_Range as Figure11_Omega_Range_string
from header import *
import pickle

try:
    from ExtraSpecifications import *
except:
    pass

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

    fig = plt.figure(figsize=(4, 4))
    plt.subplot(2,1,1)
    fig.subplots_adjust(left=0.15)
    axes = fig.gca()
    sns_map = sns.heatmap(
        data[::, ::-1].T,
        ax=axes,
        norm=log_norm,
        xticklabels=omega_range_string[::2],
        yticklabels=SNR_range_string[::-2],
        cbar_kws={"ticks": cbar_ticks[1:]},
    )
    sns_map.set_xlabel(r"Bandwidth ($\Omega$)")
    sns_map.set_ylabel("SNR")
    sns_map.set_title("Reconstruction Error")
    fig = sns_map.get_figure()
    axes = fig.gca()
    # axes.set_xticks([0, 1,2,3,4])
    plt.xticks(rotation = 0)
    plt.text(-0.3, 0.5,'(a)', ha='center', va='center', transform=axes.transAxes)
    axes.set_xticks([0.5, 2.5, 4.5, 6.5, 8.5, 10.5])
    axes.set_yticks([0.5,2.5,4.5,6.5,8.5])
    plt.yticks(rotation = 0)
    plt.ylim(9,0)
    fig.subplots_adjust(bottom=0.2)




    data_filename = Data_Path+"Figure11_b_VarSNR_VarShift.pkl"

    with open(data_filename, "rb") as f:  # Python 3: open(..., 'wb')
        obj = pickle.load(f, encoding="latin1")

    err_double = obj[0]
    err_single = obj[1]
    Figure11_noise_range = obj[2]
    SNR_range_string = obj[3]
    shifts = obj[4]
    shifts_range_string = obj[5]

    shifts_range_string.append("1-channel")
    shifts_range_string.append("1-channel")

    data = np.mean(err_double, 2)
    data_aug = np.zeros((data.shape[0], data.shape[1] + 2))
    data_aug[:, : data.shape[1]] = data[:, :]
    data_aug[:, -2] = np.mean(err_single, 1)
    data_aug[:, -1] = np.mean(err_single, 1)
    min_val = data.min().min()
    max_val = data.max().max()
    min_val = 1e-8
    max_val = 1e-1
    log_norm = LogNorm(vmin=min_val, vmax=max_val)

    cbar_ticks = [
        math.pow(10, i)
        for i in range(
            math.floor(math.log10(min_val)),
            1 + math.ceil(math.log10(max_val)),
        )
    ]

    plt.subplot(2,1,2)
    fig.subplots_adjust(left=0.15)
    axes = fig.gca()
    sns_map = sns.heatmap(
        data_aug[::-1, :],
        ax=axes,
        norm=log_norm,
        yticklabels=SNR_range_string[::-2],
        xticklabels=shifts_range_string[::2],
        cbar_kws={"ticks": cbar_ticks[1:]},
    )
    sns_map.set_ylabel(r"SNR")
    sns_map.set_xlabel("Shift")
    sns_map.set_title("Reconstruction Error")
    fig = sns_map.get_figure()
    axes = fig.gca()
    axes.set_xticks([0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18])
    # plt.yticks(rotation = 0)
    # axes.set_yticks([0, 2, 4, 6, 8, 10])
    axes.set_yticks([0.5,2.5,4.5,6.5,8.5])
    plt.yticks(rotation = 0)
    plt.ylim(9,0)
    plt.text(-0.3, 0.5,'(b)', ha='center', va='center', transform=axes.transAxes)
    axes.vlines(17, *axes.get_ylim(), color="yellow", linestyle="--")
    fig.subplots_adjust(bottom=0.2)

    plt.tight_layout()


    if To_Svg:
        figure_filename = Figure_Path+"Figure11_combined.svg"
        fig.savefig(figure_filename)
    else:
        figure_filename = Figure_Path+"Figure11_combined.png"
        fig.savefig(figure_filename, dpi=600)


if __name__ == "__main__":
    Generate()
