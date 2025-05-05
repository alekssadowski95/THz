import numpy as np
import matplotlib.pyplot as plt

# Definition of functions

# Calculates the THz beam size considering that the initial beam waist is w0 and that the propagation distance is z
def calculate_thz_beam_size(w0, z, lambda_thz):
    k_THz = 2 * np.pi / lambda_thz
    if k_THz * w0 / 2 > 5:
        z_diff = (np.pi * w0**2) / lambda_thz
    else:
        z_diff = (2 * np.pi**2 * w0**3) / lambda_thz**2
    return w0 * np.sqrt(1 + (z / z_diff)**2)

# Calculates the collection efficiency of the off-axis parabolic mirror (OAPM) knowing the beam size
def calculate_collection_efficiency(w_z, r_mirror):
    return 1 - np.exp(-2 * r_mirror**2 / w_z**2)

# Calculates the focal length of a lens used to focalize a collimated beam of width w_l to a widht of w_0 at the focus
def calculate_focal_length(w_l, w_0, lambda_IR):
    return (np.pi * w_l * w_0) / lambda_IR

# Optimize the value of w0 to maximize the collection efficiency and take the optimized w0 that has the highest THz field strength
def optimize_w_0(lambda_thz_list, lambda_IR, z_mirror, r_mirror, w_l, efficiency_threshold, w0_i, w0_f, num_points):
    wl_values = np.linspace(w0_i, w0_f, num_points)
    w0_values = wl_values / np.sqrt(2)

    plt.figure()
    for lambda_thz in lambda_thz_list:
        efficiencies = []
        w0_above_threshold = []
        for w0 in w0_values:
            wz = calculate_thz_beam_size(w0, z_mirror, lambda_thz)
            eff = calculate_collection_efficiency(wz, r_mirror) * 100
            efficiencies.append(eff)
            if eff > efficiency_threshold:
                w0_above_threshold.append(w0)
        if w0_above_threshold:
            best_w0 = w0_above_threshold[0]
            focal_length = calculate_focal_length(w_l, best_w0, lambda_IR)
            print(f"λ = {lambda_thz*1e6:.0f} µm → Optimal w0: {best_w0*1e3:.3f} mm, Focal length: {focal_length*1e3:.2f} mm")
        else:
            print(f"λ = {lambda_thz*1e6:.0f} µm → No w0 found above {efficiency_threshold}% efficiency threshold")

        plt.semilogx(w0_values * 1e3, efficiencies, label=f'λ = {lambda_thz*1e6:.0f} µm')

    plt.xlabel("Initial Beam Waist w0 (mm)")
    plt.ylabel("Collection Efficiency (%)")
    plt.grid(True)
    plt.ylim(0, 105)
    plt.axhline(y=efficiency_threshold, color='black', linestyle='--', label=f'Threshold = {efficiency_threshold}%')
    plt.legend()
    plt.show()

# Definitions of parameters
lambda_thz_list = [60e-6, 75e-6, 100e-6, 150e-6, 300e-6]    # List of tested 'THz' wavelengths (m)
lambda_IR = 800e-9                                          # Wavelength of the IR laser (m)
z_mirror = 15.0e-3                                          # Distance from the center of the generation crystal to the OAPM (m)
r_mirror = 6.35e-3                                          # Radius of the OAPM (m)
w_l = 11.0e-3                                               # Width of the focused IR laser beam (m)
efficiency_threshold = 95                                   # Minimum efficiency threshold
w0_start = 1e-6                                             # Start of the tested beam waists (m)
w0_end = 10.0e-3                                            # End of the tested beam waists (m)
nb_of_points = 10000                                        # Number of points used for the calculation

# Utilization of functions
optimize_w_0(lambda_thz_list, lambda_IR, z_mirror, r_mirror, w_l, efficiency_threshold, w0_start, w0_end, nb_of_points)