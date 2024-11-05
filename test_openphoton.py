import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import openphoton as op
import numpy as np

# Définition des fonctions
def beam_init():
    return op.devices.laser_beam(side_length=L, aperture_radius=0.05, pixels=N)

def E_0_calc(D_1, D_2, f_length, theta):
    u_0 = beam_init()
    u_1 = op.rayleigh_sommerfeld.fresnel_approx(u_0, L, wavelength, D_1/np.cos(theta))
    u_2 = np.multiply(u_1, op.lenses.converging_lens(u_1,L,wavelength,f_length))
    return op.rayleigh_sommerfeld.fresnel_approx(u_2, L, wavelength, D_2/np.cos(theta))

def E_1_calc(D_1, D_2, f_length, d, theta):
    u_0 = beam_init()
    d_r = d/np.cos(theta)
    u_1 = op.rayleigh_sommerfeld.fresnel_approx(u_0, L, wavelength, D_1)
    u_2 = np.multiply(u_1, op.lenses.converging_lens(u_1,L,wavelength,f_length))
    u_3 = op.rayleigh_sommerfeld.fresnel_approx(u_2, L, wavelength, 2 * d_r)
    u_4 = np.multiply(u_3, op.lenses.converging_lens(u_3,L,wavelength,-f_length))
    u_5 = op.rayleigh_sommerfeld.fresnel_approx(u_4, L, wavelength, 2 * d_r)
    return op.rayleigh_sommerfeld.fresnel_approx(u_5, L, wavelength, D_2)

def plot_int(E):
    intensity = np.abs(E)**2
    plt.imshow(intensity)
    plt.colorbar(label='Intensité (unité)')
    plt.title('Intensité sur le plan image')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

def plot_phase(E):
    phase = np.angle(E)
    plt.imshow(phase)
    plt.colorbar(label='Phase (rad)')
    plt.title('Phase sur le plan image')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

# Définition des constantes
L = 5e-3
N = 1024
dx = L/N
wavelength = 800e-9 #Longueur d'onde (m)

# Exécution des fonctions
E_out = E_0_calc(1,1,0.25,0) + E_1_calc(1,1,0.25,5e-3,0)
plot_int(E_out)
plot_phase(E_out)