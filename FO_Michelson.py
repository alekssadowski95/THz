# Importation of relevent libraries
import imageio.v2 as imageio
import numpy as np
import matplotlib.pyplot as plt

# Definition of functions
# Calculation of the traveled distance (laser to screen) for both paths in the Michelson interferometer
def Calc_geo(theta, D, C, abs_G, d):
    # Calculation of the relevent positions in the interferometer
    x_1 = D/(1 - np.tan(theta))
    y_1 = D * np.tan(theta)/(1 - np.tan(theta))
    x_2 = (C - y_1) * np.tan(theta) + x_1
    y_2 = C
    x_3 = (x_2 * np.tan(theta + np.pi/2) - y_2 - D)/(np.tan(theta + np.pi/2) - 1)
    y_3 = x_3 - D
    x_4 = D + C + d
    y_4 = (C + d) * np.tan(theta) + y_1
    x_5 = (x_4 * np.tan(np.pi - theta) - y_4 - D)/(np.tan(np.pi - theta) - 1)
    y_5 = x_5 - D

    # Calculation of the distances between each position
    Pl_1 = x_1/np.cos(theta)
    Pl_2_1 = (C - y_1)/np.cos(theta)
    Pl_2_2 = (C + d)/np.cos(theta)
    Pl_3_1 = (C - y_3)/np.cos(theta)
    Pl_3_2 = (x_4 - x_5)/np.cos(theta)
    Pl_4_1 = (y_3 + abs_G)/np.cos(theta)
    Pl_4_2 = (y_5 + abs_G)/np.cos(theta)

    # Calculation of the length of both optical paths
    path_1 = Pl_1 + Pl_2_1 + Pl_3_1 + Pl_4_1
    path_2 = Pl_1 + Pl_2_2 + Pl_3_2 + Pl_4_2
    return path_1, path_2

# Calculation of the output electric field associated with light traveling in one path of 
# the interferometer (Gaussian beam + free space propagation (d_n))
def E_n(d_n):
    W_d_n = W_0 * np.sqrt(1 + (d_n / z_R)**2)
    R_d_n = d_n * (1 + (z_R / d_n)**2)
    gouy_phase = np.arctan(d_n/z_R)
    tot_phase = np.exp(1j * (-k * d_n - k * (X**2 + Y**2)/(2 * R_d_n) + gouy_phase))
    envelope = np.exp(- (X**2 + Y**2) / (W_d_n**2))
    E_n_field = (A_0 * W_0 / W_d_n) * envelope * tot_phase
    return E_n_field

# Calculation of the interference pattern of the interferometer
def Michelson(theta, D, C, abs_G, d):
    [d1,d2] = Calc_geo(theta, D, C, abs_G, d)
    E_out = E_n(d1) + E_n(d2)
    return E_out

# Plotting normalized intensity
def plot_int(E, show_cond):
    intensity = np.abs(E)**2
    norm_int = intensity / np.max(intensity)
    #plt.figure(figsize=(8,6), dpi=100)
    plt.imshow(norm_int, extent=[-L_x/2, L_x/2, -L_y/2, L_y/2], 
               origin='lower', cmap='inferno', aspect='auto')
    plt.colorbar(label='Normalized Intensity')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    if show_cond == True:
        plt.show()

# Plotting phase
def plot_phase(E, show_cond):
    phase = np.angle(E)
    plt.imshow(phase, extent=[-L_x/2, L_x/2, -L_y/2, L_y/2], 
               origin='lower', cmap='inferno', aspect='auto')
    plt.colorbar(label='Phase (rad)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    if show_cond == True:
        plt.show()

# Plotting intensity profile along an axis
def plot_int_profile(E, Axis, position):
    plt.figure(figsize=(8, 6))

    intensity = np.abs(E)**2
    norm_int = intensity / np.max(intensity)
    if Axis == 'x':
        index = np.argmin(np.abs(y - position))
        I_profile = norm_int[index, :]
        plt.title("Intensity profile along the x-axis (y = {} m)".format(y[index]))
        plt.xlabel('x (m)')
    if Axis == 'y':
        index = np.argmin(np.abs(x - position))
        I_profile = norm_int[:, index]
        plt.title("Intensity profile along the y-axis (x = {} m)".format(x[index]))
        plt.xlabel('y (m)')

    plt.plot(x, I_profile)
    plt.ylabel('Normalized intensity')
    plt.show()

# Creating a GIF of the variation of the interference pattern with respect to the incident angle
def Pattern_vs_angle_GIF(init_angle, fin_angle, nb_frames):
    # Definition of the angle range
    angles = np.linspace(init_angle, fin_angle, nb_frames)

    # Initialization of the frames of the GIF
    frames = []

    # Calculation and save of the frames of the GIF
    for i in range(len(angles)):
        E_out = Michelson(angles[i], 0.2, 0.1, 0.2, 0.5)
        plot_int(E_out, False)
        plt.text(0.05, 0.95, f"Angle: {angles[i]*1e3:.4f} mrad", 
                 transform=plt.gca().transAxes, 
                 fontsize=12, color='red', 
                 verticalalignment='top', 
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
        filename = f'frame_{i+1}.png'
        plt.savefig(filename)
        frames.append(imageio.imread(filename))
        plt.close()
    
    # Combination of the frames to create the GIF
    imageio.mimsave('animation.gif', frames, fps = 5)
    print('GIF created successfully!')

# Constants
wl = 632.8e-9               # Wavelength (m)
A_0 = 1                     # Amplitude (Amplitude units)
k = 2 * np.pi / wl          # Wavenumber (m^-1)
W_0 = 1e-4                  # Minimum beam waist (m)
z_R = np.pi * W_0**2 / wl   # Rayleigh range (m)

# Grid setup
L_x = 5e-3                  # Horizontal length of the grid (m)
L_y = 5e-3                  # Vertical length of the grid (m)
grid_size = (1000, 1000)    # Grid resolution (pixels)
x = np.linspace(-L_x/2, L_x/2, grid_size[0])
y = np.linspace(-L_y/2, L_y/2, grid_size[1])
X, Y = np.meshgrid(x, y)

# Utilization of the defined functions:
E_out = Michelson(0, 0.2, 0.1, 0.2, 0.1)
plot_int(E_out, True)
# plot_int_profile(E_out, 'x', 1e-3)
# plot_int_profile(E_out, 'y', -2e-3)
# Pattern_vs_angle_GIF(-0.0002, 0.0002, 50)