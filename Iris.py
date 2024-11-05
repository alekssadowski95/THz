# Importation des librairies nécessaires
import numpy as np
import matplotlib.pyplot as plt

# Définition des fonctions
def allignement_error(d_t, d_I_1, d_I_2):
    theta_e = np.atan((2 * delta_p_y + r_I_1 + r_I_2)/(d_I_2 - d_I_1 - 2 * delta_p_x))
    delta_y = d_t * np.tan(theta_e)
    return theta_e, delta_y

def data_all_error(d_t_range, d_I_1_range, d_I_2_range):
    d_t_def = 1
    d_I_1_def = 0.4
    d_I_2_def = 0.8

    t_te = np.zeros(len(d_t_range))
    t_dy = np.zeros(len(d_t_range))
    I_1_te = np.zeros(len(d_I_1_range))
    I_1_dy = np.zeros(len(d_I_1_range))
    I_2_te = np.zeros(len(d_I_2_range))
    I_2_dy = np.zeros(len(d_I_2_range))


    for i in range(len(d_t_range)):
        [t_te[i], t_dy[i]] = allignement_error(d_t_range[i], d_I_1_def, d_I_2_def)
    for o in range(len(d_I_1_range)):
        [I_1_te[o], I_1_dy[o]] = allignement_error(d_t_def, d_I_1_range[o], d_I_2_def)
    for n in range(len(d_I_2_range)):
        [I_2_te[n], I_2_dy[n]] = allignement_error(d_t_def, d_I_1_def, d_I_2_range[n])
    
    return t_te, t_dy, I_1_te, I_1_dy, I_2_te, I_2_dy

def m_set_dist_error(dist_M1_M2_range):
    te = np.zeros(len(dist_M1_M2_range))
    dy = np.zeros(len(dist_M1_M2_range))
    for i in range(len(dist_M1_M2_range)):
        d_I_1 = 0.01
        d_I_2 = d_I_1 + dist_M1_M2_range[i]
        d_t = d_I_2 + 0.01
        [te[i], dy[i]] = allignement_error(d_t, d_I_1, d_I_2)
    return te, dy

def display_error(x, y1, y2, x_lab):
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    ax1.scatter(x, y1)
    ax1.set_xlabel(x_lab)
    ax1.set_ylabel('Erreur angulaire (rad)')
    ax1.grid(True)
    ax2.scatter(x, y2)
    ax2.set_xlabel(x_lab)
    ax2.set_ylabel("Erreur linéaire sur l'axe y (m)")
    ax2.grid(True)

    plt.tight_layout()
    plt.show()

# Définition des constantes
delta_p_y = 3e-3                            # Désalignement selon l'axe y (m)
delta_p_x = 3e-3                            # Désalignement selon l'axe x (m)
r_I_1 = 2e-3                                # Rayon de l'iris 1
r_I_2 = 2e-3                                # Rayon de l'iris 2
d_t_list = np.linspace(0.805, 2, 100)       # Liste de valeurs de d_t (m)
d_I_1_list = np.linspace(0.01, 0.795, 100)  # Liste de valeurs de d_t (m)
d_I_2_list = np.linspace(0.405, 0.995, 100) # Liste de valeurs de d_t (m)
m_dist_list = np.linspace(0.1, 1, 100)      # Liste de distances entre les 2 miroirs (m)

error_data = data_all_error(d_t_list, d_I_1_list, d_I_2_list)
err_fix_m_pos = m_set_dist_error(m_dist_list)
#display_error(d_t_list, error_data[0], error_data[1], 'Distance totale (m)')
#display_error(d_I_1_list, error_data[2], error_data[3], 'Distance origine-miroir 1 (m)')
#display_error(d_I_2_list, error_data[4], error_data[5], 'Distance origine-miroir 2 (m)')
display_error(m_dist_list, err_fix_m_pos[0], err_fix_m_pos[1], 'Distance miroir 1-miroir 2 (m)')
