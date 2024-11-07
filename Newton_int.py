# Importation of libraries
import numpy as np
import imageio.v2 as imageio
import matplotlib.pyplot as plt

# Definition of functions

# Calculation of t (Added distance between the two glass components at radial distance r from the center 
# of the convexe lens due to its curvature)
def calc_t(r):
    return r**2 / (2 * R)

# Calculation of the output intensity of the optical system
def calc_int(theta):
    r = np.sqrt(X**2 + Y**2)
    phase_difference = 2 * np.pi / wl * (2 * np.cos(theta) * (calc_t(r) + tau) + 2 * epsilon)
    intensity = I0 / (1 + np.sin(phase_difference/2)**2)
    return intensity

# Usage of calc_int to display the output intensity of the optical system
def plot_int_with_smoothing(int, theta, sigma, cond_show_plt):
    plt.figure(figsize=(8, 8))
    from scipy.ndimage import gaussian_filter
    smoothed_intensity = gaussian_filter(int, sigma)
    plt.imshow(smoothed_intensity, cmap='inferno', interpolation='bilinear')
    plt.colorbar(label="Intensity")
    plt.title(f"Newton Interference Pattern\nAngle of Incidence = {np.degrees(theta):.1f}°")
    plt.xlabel("X (meters)")
    plt.ylabel("Y (meters)")
    plt.axis('equal')
    if cond_show_plt == True:
        plt.show()

# Creates an animation with the output of plot_int_with_smoothing for different values of incident angles
def animate_res_angles(init_angle, fin_angle, nb_frames):
    # Définitions des angles qui seront évalué
    angles = np.linspace(init_angle, fin_angle, nb_frames)

    # Initialisation des frames du gif
    frames = []

    # Calculs et sauvegarde des images du gif
    for i in range(len(angles)):
        plot_int_with_smoothing(calc_int(angles[i]), angles[i], 1, False)
        filename = f'frame_{i+1}.png'
        plt.savefig(filename)
        frames.append(imageio.imread(filename))
        plt.close()
    
    # Création d'un .gif à partir des images sauvegardées
    imageio.mimsave('animation.gif', frames, duration = 5)
    print('GIF animation created successfully!')

# Constants
wl = 632.8e-9       # Wavelength (m)
R = 103e-3          # Radius of curvature of the lens (m)
I0 = 1              # Maximum intensity (a.u.)
d_i_lens = 1        # Distance between the laser and the lens (m)
tau = 0             # Distance between the two glass components (m)
epsilon = 0         # Effective path change of a single reflection (m)
screen_size = 4e-3  # Size of screen (meters)
num_points = 1000   # Number of grid points per axis

# Grid setup
x = np.linspace(-screen_size/2, screen_size/2, num_points)
y = np.linspace(-screen_size/2, screen_size/2, num_points)
X, Y = np.meshgrid(x, y)

# Execution of functions
animate_res_angles(0, np.pi/4, 100)