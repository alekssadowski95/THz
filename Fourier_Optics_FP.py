# Importation des librairies nécessaires
import numpy as np
import scipy.fft as fft
import scipy.signal as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad_vec
from scipy.special import jv

# Définition des fonctions
# E(x, y, z) d'un faisceau gaussien 
def initial_signal(x, y, z):
    A_0 = 1
    z_0 = np.pi * W_0**2/wl
    W = W_0 * np.sqrt(1 + (z/z_0)**2)
    eta = np.atan(z/z_0)
    if z == 0:
        R = np.inf
    else:
        R = z * (1 + (z_0/z)**2)
    U = A_0 * W_0/W * np.exp(-(x**2 + y**2)/W**2) * np.exp(-1j*(k*z + k*(x**2 + y**2)/(2*R) - eta))
    return U

# Réponse impulsionelle d'une lentille + propagation libre
def h_lens(x, y, d_1, d_2, f):

    def integrand(r, rho_val, d_1, d_2, f):
        fc_err = 1/d_1 + 1/d_2 - 1/f
        exponent = -1j * np.pi * fc_err * r**2 / wl
        bessel_term = 2 * np.pi * jv(0, 2 * np.pi * rho_val * r)
        return r * np.exp(exponent) * bessel_term
    
    h_1 = 1j / (wl * d_1) * np.exp(-1j * k * d_1)
    h_2 = 1j / (wl * d_2) * np.exp(-1j * k * d_2)
    phase_factor = np.exp(-1j * np.pi * (x**2 + y**2) / (wl * d_2))
    rho = np.sqrt(x**2 + y**2)
    P1_values, error = quad_vec(lambda r: integrand(r, rho, d_1, d_2, f), 0, ra, limit=100)
    modified_P1_values = P1_values * h_1 * h_2 * phase_factor
    
    return modified_P1_values

# Utilisation des fonctions précédentes pour calculer E(x,y,z) au plan image
def final_signals():
    # Calcul des expressions à convolué
    E_in = initial_signal(X, Y, d_i_to_M1/np.cos(theta) - d_shift)
    h_l1_0 = h_lens(X, Y, d_shift, (d + d_M2_to_f)/np.cos(theta), focal)
    h_l1_1 = h_lens(X, Y, d_shift, d/np.cos(theta), focal)
    h_l2 = h_lens(X, Y, d/np.cos(theta), (d + d_M2_to_f)/np.cos(theta), -focal)
    
    # Calcul des signaux finaux (convolutions)
    E_i_to_M = sp.convolve2d(E_in, h_l1_1, mode = 'same')
    E_0 = sp.convolve2d(E_in, h_l1_0, mode = 'same')
    E_1 = sp.convolve2d(E_i_to_M, h_l2, mode = 'same')
    E_out = E_0 + E_1

    return E_0, E_1, E_out

# Utilisation de la fonction final_signals pour afficher l'intensité du résultat
def plot_int(E):
    intensity = np.abs(E)**2
    plt.contourf(x, y, intensity, 100)
    plt.colorbar(label='Intensité (unité)')
    plt.title('Intensité sur le plan image')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

# Utilisation de la fonction final_signals pour afficher la phase du résultat
def plot_phase(E):
    phase = np.angle(E)
    plt.contourf(X, Y, phase, 100)
    plt.colorbar(label='Phase (rad)')
    plt.title('Phase sur le plan image')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

# Définition des constantes
theta = 0               # Angle d'incidence (rad)
wl = 800e-9             # Longueur d'onde (m)
k = 2 * np.pi/wl        # Nombre d'onde (m^-1)
W_0 = 10e-6             # Minimum beam waist (m)
d_i_to_M1 = 0.5         # Distance entre l'origine du faisceau et le premier miroir (m)
d_M2_to_f = 0.5         # Distance entre le deuxième miroir et l'écran (m)
ra = 12.7e-3            # Rayon de la lentille (miroir courbé)
focal = 0.025           # Distance focale du miroir (m)
d = 6.4e-3              # Distance entre les 2 miroirs (m)
d_shift = 0.05          # Paramètre de décalage de la position de calcul de E_in et de h_l1 pour éviter une division par 0
L = 5e-3                # Longueur d'un côté du détecteur (m)
grid_size = (100,100)   # Nombres de points dans le mesh

# Definition de la plage x,y
x = np.linspace(-L, L, grid_size[0])
y = np.linspace(-L, L, grid_size[1])
X, Y = np.meshgrid(x, y)

output_signals = final_signals()
plot_int(output_signals[2])

# Problème : Différents d_shift donnent différentes solutions (ne devrait pas être le cas) 
# - Tests différent d_shift:
p1 = False
if p1 == True:
    d_shifts = [0.025, 0.05, 0.1]
    for n in range(len(d_shifts)):
        d_shift = d_shifts[n]
        gauss = initial_signal(X, Y, d_i_to_M1 - d_shift)
        pup = h_lens(X, Y, d_shift, d/np.cos(theta), focal)
        conv = sp.convolve2d(gauss, pup, 'same')
        plot_int(conv)


# Fonctions utilisées précédements:
# Impulse response de la propagation libre (sous l'approximation de Fresnel)
def h_free_space_prop(x, y, p_dist):
    h_0 = 1j/(wl * p_dist) * np.exp(-1j * k * p_dist)
    h = h_0 * np.exp(-1j * k * (x**2 + y**2)/(2 * p_dist))
    return h

# TF de la fonction pupille généralisé évalué à (x/(lambda * d2), y/(lambda * d2))
def FT_gen_pupil_fcn(x, y, d_1, d_2, f):
    fc_err = 1/d_1 + 1/d_2 - 1/f
    p = (x**2 + y**2 <= ra**2) * np.exp(-1j * np.pi * fc_err * (x**2 + y**2) / wl)
    return fft.fftshift(fft.fft2(p))

# Impulse response de la propagation à travers une lentielle
def h_lens_v1(x, y, d_1, d_2, f):
    h_1 = 1j/(wl * d_1) * np.exp(-1j * k * d_1)
    h_2 = 1j/(wl * d_2) * np.exp(-1j * k * d_2)
    pup = FT_gen_pupil_fcn(x, y, d_1, d_2, f)
    h = h_1 * h_2 * np.exp(-1j * np.pi * (x**2 + y**2)/(wl * d_2)) * pup
    return h

# Utilisation des fonctions précédentes pour calculer E(x,y,z) au plan image
def final_signals_v1():
    # Calcul des expressions à convolué
    E_in = initial_signal(X, Y, d_i_to_M1/np.cos(theta) - d_shift)
    h_p =  h_free_space_prop(X, Y, d_M2_to_f/np.cos(theta))
    h_l1 = h_lens(X, Y, d_shift, d/np.cos(theta), focal)
    h_l2 = h_lens(X, Y, d/np.cos(theta), d/np.cos(theta), -focal)
    
    # Calcul des signaux finaux (convolutions)
    E_1st_lens = sp.convolve2d(E_in, h_l1, mode = 'same')
    E_2nd_lens = sp.convolve2d(E_1st_lens, h_l2, mode = 'same')
    E_0 = sp.convolve2d(E_1st_lens, h_p, mode = 'same')
    E_1 = sp.convolve2d(E_2nd_lens, h_p, mode = 'same')
    E_out = E_0 + E_1

    return E_0, E_1, E_out
