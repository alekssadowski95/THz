# Importation of relevent libraries
import imageio.v2 as imageio
import numpy as np
import matplotlib.pyplot as plt

# Definition of functions
# Calculation of the traveled distance (laser to screen) for both paths in the Michelson interferometer
def Calc_geo(theta, D, C, G, d):
    # Calculation of the relevent positions in the interferometer
    theta_pp = -(np.pi/2 - theta - M1_tilt)
    x_1 = (D + Beam_y_err - Beam_x_err * np.tan(theta)) / (1 - np.tan(theta))
    y_1 = x_1 - D
    x_2 = (y_1 + (D + M1_x_err) * np.tan(M1_tilt) - (C + M1_y_err) - x_1 * np.tan(np.pi/2 - theta)) / (np.tan(M1_tilt) - np.tan(np.pi/2 - theta))
    y_2 = (x_2 - D) * np.tan(M1_tilt) + C + M1_y_err
    x_3 = (D + y_2 - x_2 * np.tan(theta_pp)) / (1 - np.tan(theta_pp))
    y_3 = x_3 - D
    x_4 = ((D + C + d + M2_x_err) * np.tan(np.pi/2 - M2_tilt) + Beam_y_err - Beam_x_err * np.tan(theta) + M2_y_err) / (np.tan(np.pi/2 - M2_tilt) - np.tan(theta))
    y_4 = (x_4 - Beam_x_err) * np.tan(theta) + Beam_y_err
    x_5 = (y_4 + D + x_4 * np.tan(theta + M2_tilt)) / (1 + np.tan(theta + M2_tilt))
    y_5 = x_5 - D
    x_fin_1 = ((D + Det_x_err) * np.tan(Det_tilt) + G + y_2 - Det_y_err - x_2 * np.tan(theta_pp)) / (np.tan(Det_tilt) - np.tan(theta_pp))
    y_fin_1 = (x_fin_1 - (D + Det_x_err)) * np.tan(Det_tilt) - G + Det_y_err
    x_fin_2 = ((D + Det_x_err) * np.tan(Det_tilt) + G + y_5 - Det_y_err - x_5 * np.tan(np.pi/2 - theta - M2_tilt)) / (np.tan(Det_tilt) - np.tan(np.pi/2 - theta - M2_tilt))
    y_fin_2 = (x_fin_2 - (D + Det_x_err)) * np.tan(Det_tilt) - G + Det_y_err
    
    positions = [[x_1, y_1], [x_2, y_2], [x_3, y_3], [x_4, y_4], [x_5, y_5], [x_fin_1, y_fin_1], [x_fin_2, y_fin_2]]

    # Calculation of the distances between each position
    d_1 = np.sqrt(x_1**2 + y_1**2)
    d_2 = np.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2)
    d_3 = np.sqrt((x_3 - x_2)**2 + (y_3 - y_2)**2)
    d_fin_1 = np.sqrt((x_fin_1 - x_3)**2 + (y_fin_1 - y_3)**2)
    d_4 = np.sqrt((x_4 - x_1)**2 + (y_4 - y_1)**2)
    d_5 = np.sqrt((x_5 - x_4)**2 + (y_5 - y_4)**2)
    d_fin_2 = np.sqrt((x_fin_2 - x_5)**2 + (y_fin_2 - y_5)**2)

    # Calculation of the length of both optical paths
    path_1 = d_1 + d_2 + d_3 + d_fin_1
    path_2 = d_1 + d_4 + d_5 + d_fin_2
    return path_1, path_2, positions

# Displays a cartoon representation of the otpical system
def plot_geo(theta, D, C, abs_G, d):
    [_,_,pos_2] = Calc_geo(theta, D, C, abs_G, d)
    x_BS = np.linspace(D - D/4, D + D/4, 10)
    y_BS = x_BS - D
    x_Det = np.linspace((D - D/4) + Det_x_err, (D + D/4) + Det_x_err, 10)
    y_Det = (x_Det - D) * np.tan(Det_tilt) - abs_G + Det_y_err
    x_M1 = np.linspace((D - D/4) + M1_x_err, (D + D/4) + M1_x_err, 10)
    y_M1 = (x_M1 - D) * np.tan(M1_tilt) + C + M1_y_err
    y_M2 = np.linspace(-D/4 + M2_y_err, D/4 + M2_y_err, 10)
    x_M2 = y_M2/np.tan(np.pi/2 - M2_tilt) + (D + C + d) + M2_x_err

    for i in range(len(pos_2)):
        plt.scatter(pos_2[i][0], pos_2[i][1], color = 'Black')
    plt.plot([Beam_x_err, pos_2[0][0]], [Beam_y_err, pos_2[0][1]], color = 'Blue')
    plt.plot([pos_2[0][0], pos_2[1][0]], [pos_2[0][1], pos_2[1][1]], color = 'Blue')
    plt.plot([pos_2[1][0], pos_2[2][0]], [pos_2[1][1], pos_2[2][1]], color = 'Blue')
    plt.plot([pos_2[2][0], pos_2[5][0]], [pos_2[2][1], pos_2[5][1]], color = 'Blue')

    plt.plot([Beam_x_err, pos_2[0][0]], [Beam_y_err, pos_2[0][1]], color = 'Blue')
    plt.plot([pos_2[0][0], pos_2[3][0]], [pos_2[0][1], pos_2[3][1]], color = 'red')
    plt.plot([pos_2[3][0], pos_2[4][0]], [pos_2[3][1], pos_2[4][1]], color = 'red')
    plt.plot([pos_2[4][0], pos_2[6][0]], [pos_2[4][1], pos_2[6][1]], color = 'red')

    plt.plot(x_BS, y_BS, color = 'Black', label = 'Beam splitter', linestyle = 'dashed')
    plt.plot(x_Det, y_Det, color = 'Black', label = 'Detctor')
    plt.plot(x_M1, y_M1, color = 'Grey', label = 'Mirror 1')
    plt.plot(x_M2, y_M2, color = 'Green', label = 'Mirror 2')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.show()

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
    [d1,d2,_] = Calc_geo(theta, D, C, abs_G, d)
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
    imageio.mimsave('animation.gif', frames, fps = 20)
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

# Position error of optical components
M1_tilt = 0
M2_tilt = 0
Det_tilt = 0
M1_x_err = 0
M1_y_err = 0
M2_x_err = 0
M2_y_err = 0
Det_x_err = 0
Det_y_err = 0
Beam_x_err = 0
Beam_y_err = 0

# Utilization of the defined functions:
plot_geo(0.02, 0.2, 0.1, 0.2, 0.1)

E_out = Michelson(0.02, 0.2, 0.1, 0.2, 0.1)
plot_int(E_out, True)

# plot_int_profile(E_out, 'x', 1e-3)
# plot_int_profile(E_out, 'y', -2e-3)

# Pattern_vs_angle_GIF(-5e-6, 5e-6, 100)