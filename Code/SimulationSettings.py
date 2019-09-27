import numpy as np

graphical_import = True

# should be 100
n_trials = 100
n_omega = 12
#should be 10
n_shifts = 20
# n_shifts = 4


omega_range = [0.25 * np.pi]
omega_range_string = [
    r"$0.25\pi$",
    r"$0.25\sqrt{2}\pi$",
    r"$0.5\pi$",
    r"$0.5\sqrt{2}\pi$",
    r"$\pi$",
    r"$\sqrt{2}\pi$",
    r"$2\pi$",
    r"$2\sqrt{2}\pi$",
    r"$4\pi$",
    r"$4\sqrt{2}\pi$",
    r"$8\pi$",
    r"$8\sqrt{2}\pi$",
    r"$16\pi$",
]
for i in range(n_omega):
    omega_range.append(omega_range[-1] * np.sqrt(2))

shifts = [1]
shifts_range_string = [r"$1$"]
for i in range(n_shifts):
    shifts.append(shifts[-1] / np.sqrt(10))
    if i == 0:
        shifts_range_string.append(r"$\sqrt{{10}}$")
    elif i % 2 == 1:
        shifts_range_string.append(r"$10^{0}$".format(-int((i+1) / 2)))
    else:
        shifts_range_string.append(r"$10^{0}\sqrt{{10}}$".format(-int(i / 2)))

num_channels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
omega_range_M_channels = np.arange(1 * np.pi, 21 * np.pi, np.pi)
# omega_range_M_channels = np.arange(1 * np.pi, 5 * np.pi, np.pi)
omega_range_M_channels_string = [
    r"$\pi$",
    r"$2\pi$",
    r"$3\pi$",
    r"$4\pi$",
    r"$5\pi$",
    r"$6\pi$",
    r"$7\pi$",
    r"$8\pi$",
    r"$9\pi$",
    r"$10\pi$",
    r"$11\pi$",
    r"$12\pi$",
    r"$13\pi$",
    r"$14\pi$",
    r"$15\pi$",
    r"$16\pi$",
    r"$17\pi$",
    r"$18\pi$",
    r"$19\pi$",
    r"$20\pi$",
]

delta_t = 1e-4
end_time = 20