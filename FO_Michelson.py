# Importation des librairies nécessaires
import numpy as np
import scipy.fft as fft
import imageio.v2 as imageio
import matplotlib.pyplot as plt

# Définition des fonctions
# Transfromé de Fourier d'un faisceau gaussien
def Gaussian_beam_FD():
    E = A_0 * np.exp(-(X**2 + Y**2) / (2 * W_0**2))
    E_FD = fft.fftshift(fft.fft2(E))
    return E_FD

# Calcul de la fonction de transfert de l'espace libre (free space propagation)
def Transfer_func_FSP(d):
    args = np.maximum(1 / wl**2 - (X/(wl * d))**2 - (Y/(wl * d))**2, 0)
    H = np.exp(-1j * 2 * np.pi * np.sqrt(args) * d)
    return H

# Calcul de la distance parcourue par la lumière dans chaque partie de l'interféromètre 
def Calc_geo(theta, D, C, abs_G, d):
    # Calculs des positions
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

    # Calculs des distances entre chaque réflexions/transmission
    Pl_1 = x_1/np.cos(theta)
    Pl_2_1 = (C - y_1)/np.cos(theta)
    Pl_2_2 = (C + d)/np.cos(theta)
    Pl_3_1 = (C - y_3)/np.cos(theta)
    Pl_3_2 = (x_4 - x_5)/np.cos(theta)
    Pl_4_1 = (y_3 + abs_G)/np.cos(theta)
    Pl_4_2 = (y_5 + abs_G)/np.cos(theta)

    # Calculs de la longueur de chaque parcours optique
    path_1 = Pl_1 + Pl_2_1 + Pl_3_1 + Pl_4_1
    path_2 = Pl_1 + Pl_2_2 + Pl_3_2 + Pl_4_2
    return path_1, path_2

# Calcul du champ électrique final
def E_calc(theta, D, C, abs_G, d):
    # Calculs des distances pour la configuration du système donné en entré
    d_n = Calc_geo(theta, D, C, abs_G, d)

    # Calculs des fonctions de transferts individuelles
    H_in = Gaussian_beam_FD()
    H_1 = Transfer_func_FSP(d_n[0])
    H_2 = Transfer_func_FSP(d_n[1])

    # Calculs des fonctions de transferts de chaque chemin
    FD_E_1 = H_in * H_1
    FD_E_2 = H_in * H_2
    
    # Calculs du champ électrique à partir des fonctions de transferts (FFT)
    E_out = np.fft.ifft2(fft.ifftshift(FD_E_1 + FD_E_2))
    return E_out

# Affichage de l'intensité du champ électrique final
def plot_int(E):
    intensity = np.abs(E)**2
    plt.contourf(X, Y, intensity, 500, cmap = 'inferno')
    plt.colorbar(label='Intensité (unité)')
    plt.title('Intensité sur le plan image')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

# Affichage de la phase du champ électrique final
def plot_phase(E):
    phase = np.angle(E)
    plt.contourf(X, Y, phase, 500, cmap ='inferno')
    plt.colorbar(label='Phase (rad)')
    plt.title('Phase sur le plan image')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

# Affichage de l'intensité du champ électrique final en fonction de x pour une valeur de y donnée
def plot_int_at_y(E, y_value):
    intensity = np.abs(E)**2
    intensity_cross_section = intensity[np.argmin(np.abs(y - y_value)), :]
    plt.plot(x, intensity_cross_section)
    plt.xlabel('x (m)')
    plt.ylabel('Intensity (add units)')
    plt.grid(True)
    plt.show()


# Définition des constantes
wl = 632.8e-9             # Longueur d'onde (m)
A_0 = 1                   # Intensité du faisceau laser (unité)
k = 2 * np.pi/wl          # Nombre d'onde (m^-1)
W_0 = 5e-4                # Minimum beam waist (m)
L_x = 5e-3                # Largeur du détecteur (m)
L_y = 5e-3                # Hauteur du détecteur (m)
grid_size = (2000,2000)   # Nombres de points dans le mesh

# Definition de la plage x,y
x = np.linspace(-L_x/2, L_x/2, grid_size[0])
y = np.linspace(-L_y/2, L_y/2, grid_size[1])
X, Y = np.meshgrid(x, y)

# Utilisation des fonctions
E_out = E_calc(0, 0.5, 0.5, 0.5, 1e-3)
plot_int(E_out)
