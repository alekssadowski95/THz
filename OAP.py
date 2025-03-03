import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

def rotate_coordinates(x, y, theta, x_center, y_center):
    x_shifted = x - x_center
    y_shifted = y - y_center

    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)], 
                                [np.sin(theta), np.cos(theta)]])
    rotated_coords = np.dot(rotation_matrix, np.array([x_shifted, y_shifted]))
    return rotated_coords[0] + x_center, rotated_coords[1] + y_center

def parent_parabola(x, f):
    return (x**2) / (4 * f) - f

def off_axis_parabolic_mirror(x, f, d_offset, invert, y_offset, x_shift, rotation_angle, flt_or_array, diff_centers_x, diff_centers_y):
    # Define the initial off axis parabolic mirror
    y = parent_parabola(x + d_offset, f) - y_offset
    if x_shift != 0:
        x = x + x_shift
    x_center = d_offset
    y_center = parent_parabola(d_offset, f) - y_offset

    # Calculate the center position of the parabola
    if flt_or_array == 'array':
        x_mid = x[int(len(x)/2)]
        y_mid = y[int(len(y)/2)]

    # Conditionnal rotation of the parabola
    if rotation_angle != 0:
        x, y = rotate_coordinates(x, y, rotation_angle, x_center, y_center)
        # shifting of the parabola so its center does not move due to the rotation
        if flt_or_array == 'array':
            diff_centers_x = x_mid - x[int(len(x)/2)]
            diff_centers_y = y_mid - y[int(len(x)/2)]
        x += diff_centers_x
        y += diff_centers_y

    
    # Conditionnal inversion of the parabola
    if invert:
        y = -y

    return x, y, diff_centers_x, diff_centers_y

def find_beam_mirror_int(start, direction, f, d_offset, invert, y_offset, x_shift, rot, diff_x, diff_y):
    min_diff = float('inf')  # Initialize the minimum difference as infinity
    best_x, best_y = None, None  # Variables to store the best x, y corresponding to the smallest difference
    for t in np.linspace(0, 100, 10000):  # Consider adjusting the range or step size if necessary
        x = start[0] + t * direction[0]
        y = start[1] + t * direction[1]
        mirror_out = off_axis_parabolic_mirror(x - x_shift, f, d_offset, invert, y_offset, x_shift, rot, 'flt', diff_x, diff_y)
        mirror_x = mirror_out[0]
        mirror_y = mirror_out[1]
        diff = np.abs(y - mirror_y)  # Calculate the absolute difference between y and mirror_y
        if diff < min_diff:  # If this difference is the smallest so far, update the best x and y
            min_diff = diff
            best_x, best_y = mirror_x, mirror_y
    return best_x, best_y  # Return the x and y corresponding to the smallest difference

def normal_vector(vals, x_query):
    x_vals = vals[0]
    y_vals = vals[1]

    # Compute numerical derivatives (tangent vector)
    dx = np.gradient(x_vals)
    dy = np.gradient(y_vals)

    # Find the closest index to x_query
    idx = np.argmin(np.abs(x_vals - x_query))

    # Compute the normal as perpendicular to the local tangent
    normal_x = -dy[idx]
    normal_y = dx[idx]

    # Normalize the normal vector
    norm = np.sqrt(normal_x**2 + normal_y**2)
    normal = np.array([normal_x / norm, normal_y / norm], dtype=np.float64)
    return normal

def reflect_ray(incident, normal):
    return incident - 2 * np.dot(incident, normal) * normal

def reflection(beam_start, beam_direction, p_focal, r_focal, inv_cond, y_offset, first_reflexion, x_shift, rot, parabola, diff_x, diff_y, Source):
    # Calculation of the intersection between the incident beam and the mirror using the direction and the origin of the incident beam
    int_x, int_y = find_beam_mirror_int(beam_start, beam_direction, p_focal, r_focal, inv_cond, y_offset, x_shift, rot, diff_x, diff_y)

    # Definition of the incident direction based on the position of the intersection
    if first_reflexion == True:
        incident_dir = np.array([int_x - Source[0], int_y - Source[1]])
        incident_dir /= np.linalg.norm(incident_dir, axis=0)
    else :
        norm_of_IncBeam = np.sqrt((int_x - beam_start[0])**2 + (int_y - beam_start[1])**2)
        incident_dir = [(int_x - beam_start[0])/norm_of_IncBeam, (int_y - beam_start[1])/norm_of_IncBeam]

    # Calculation of the normal vector at the intersection and the direction of the reflected beam
    normal_vec = normal_vector(parabola, int_x)
    r_direction = reflect_ray(incident_dir, normal_vec)

    return int_x, int_y, r_direction

def int_of_2_lines(beam_1_points, beam_2_points):
    a_1 = (beam_1_points[1][1] - beam_1_points[0][1]) / (beam_1_points[1][0] - beam_1_points[0][0])
    b_1 = beam_1_points[1][1] - a_1 * beam_1_points[1][0]
    a_2 = (beam_2_points[1][1] - beam_2_points[0][1]) / (beam_2_points[1][0] - beam_2_points[0][0])
    b_2 = beam_2_points[1][1] - a_2 * beam_2_points[1][0]

    x_int = (b_2 - b_1) / (a_1 - a_2)
    y_int = a_1 * x_int + b_1
    return x_int, y_int

def avg_distance(beams):
    beam_1 = beams[0]
    beam_2 = beams[1]
    beam_3 = beams[2]

    aligned_focal_point = [r_focal_2, v_offset_M2_0]

    intersection_1 = int_of_2_lines(beam_1, beam_2)
    intersection_2 = int_of_2_lines(beam_2, beam_3)
    intersection_3 = int_of_2_lines(beam_3, beam_1)

    dist_1 = np.sqrt(np.abs(intersection_1[0] - aligned_focal_point[0])**2 + np.abs(intersection_1[1] - aligned_focal_point[1])**2)
    dist_2 = np.sqrt(np.abs(intersection_2[0] - aligned_focal_point[0])**2 + np.abs(intersection_2[1] - aligned_focal_point[1])**2)
    dist_3 = np.sqrt(np.abs(intersection_3[0] - aligned_focal_point[0])**2 + np.abs(intersection_3[1] - aligned_focal_point[1])**2)
    avg_dist = (dist_1 + dist_2 + dist_3) / 3

    return intersection_1, intersection_2, intersection_3, avg_dist

def setup_full_calc(x1, y1, theta_1, x2, y2, theta_2, source_offset_1, source_offset_2, div, plt_cond):
    # Define shifted variables
    Source = [Source_0[0] + source_offset_1, Source_0[1] + source_offset_2]
    div_angle = div_angle_0 + div
    v_offset_M1 = v_offset_M1_0 + y1
    v_offset_M2 = v_offset_M2_0 + y2                
    h_offset_M1 = h_offset_M1_0 + x1            
    h_offset_M2 = h_offset_M2_0 + x2     
    rot_M1 = rot_M1_0 + theta_1                      
    rot_M2 = rot_M2_0 + theta_2     

    # Define the off-axis parabolic mirrors (OAP)
    parabola_1 = off_axis_parabolic_mirror(x, p_focal_1, r_focal_1, False, v_offset_M1, h_offset_M1, np.deg2rad(rot_M1), 'array', 0, 0)
    parabola_2 = off_axis_parabolic_mirror(x, p_focal_2, -r_focal_2, True, v_offset_M2, h_offset_M2, np.deg2rad(rot_M2), 'array', 0, 0)

    # Define the incident beam (The beam is illustrated by three lines : One is representing the center of the beam and the 
    # two others represent its extremities. Each line is defined with its origin and its angle with respect to the x axis).
    incident_beams = [[Source, [1, np.tan(np.deg2rad(div_angle)/2)]],
                    [Source, [1, 0]],
                    [Source, [1, -np.tan(np.deg2rad(div_angle)/2)]]]

    # Define a variable to store the start and end points of each final beam
    final_beams = [[[0, 0], [0, 0]], [[0, 0], [0, 0]], [[0, 0], [0, 0]]]  # Structure : [[Beam_1_start, Beam_1_end], [...], [...]]

    # Calculate and plot the beam path
    for i in range(len(incident_beams)):
        # Calculation of the important intersections and directions of the beam
        [int_1_x, int_1_y, dir_r1] = reflection(incident_beams[i][0], incident_beams[i][1], p_focal_1, r_focal_1, False, v_offset_M1, True, h_offset_M1, np.deg2rad(rot_M1), parabola_1, parabola_1[2], parabola_1[3], Source)
        [int_x_2, int_y_2, dir_r2] = reflection((int_1_x, int_1_y), dir_r1, p_focal_2, -r_focal_2, True, v_offset_M2, False, h_offset_M2, np.deg2rad(rot_M2), parabola_2, parabola_2[2], parabola_2[3], Source)

        # Definition of the last section of the beam path (After the second reflection)
        reflect_x_end = int_x_2 + dir_r2[0] * 40
        reflect_y_end = int_y_2 + dir_r2[1] * 40

        # Insert data in final_beams
        final_beams[i] = [[int_x_2, int_y_2], [reflect_x_end, reflect_y_end]]

        # Plot
        if plt_cond == True:
            plt.scatter(int_1_x, int_1_y, color='Blue', label= 'Intersection : Beam and 1st OAP' if i == 1 else '')
            plt.scatter(int_x_2, int_y_2, color='Grey', label='Intersection : Beam and 2nd OAP' if i == 1 else '')
            plt.plot([Source[0], int_1_x], [Source[1], int_1_y], color='Red', label='Beam' if i == 1 else '')
            plt.plot([int_1_x, int_x_2], [int_1_y, int_y_2], color='Red')
            plt.plot([int_x_2, reflect_x_end], [int_y_2, reflect_y_end], color='Red')

    crossings = avg_distance(final_beams)

    # Display the plot
    if plt_cond == True:
        plt.scatter(crossings[0][0], crossings[0][1])
        plt.scatter(crossings[1][0], crossings[1][1])
        plt.scatter(crossings[2][0], crossings[2][1])
        plt.plot(parabola_1[0], parabola_1[1], color='Black', label='OAP 1')
        plt.plot(parabola_2[0], parabola_2[1], color='Black', label='OAP 2')
        plt.scatter(25.4, 50, color='Purple', label='Focal point of the aligned system')
        plt.legend()
        plt.show()

    return crossings[3]

# Define base variables
x = np.linspace(-10, 10, 10000) # Define a range of values
Source_0 = [-15, 0]             # Position of the initial point of divergence of the beam
div_angle_0 = 25                # Angle of divergence of the initial beam in degrees (Angle between the top and the bottom of the beam)
p_focal_1 = 7.5                 # Parent focal length of the first mirror in mm
r_focal_1 = 15                  # Reflected focal length of the first mirror in mm
p_focal_2 = 12.7                # Parent focal length of the second mirror in mm
r_focal_2 = 25.4                # Reflected focal length of the second mirror in mm
v_offset_M1_0 = 0               # Vertical offset of the first mirror in mm
v_offset_M2_0 = 50              # Vertical offset of the second mirror in mm
h_offset_M1_0 = 0               # Horizontal offset of the first mirror
h_offset_M2_0 = 0               # Horizontal offset of the second mirror
rot_M1_0 = 0                    # Rotation of the first OAP in degrees
rot_M2_0 = 0                    # Rotation of the second OAP in degrees


angular_displacement = np.array([1e-4, 1e-3, 1e-2, 1e-1, 1])
lin_displacement = np.array([1e-4, 1e-3, 1e-2, 1e-1, 1])
component_variation = ['M1 horizontal displacement', 'M1 vertical displacement', 'M1 angular displacement', 
                       'M2 horizontal displacement', 'M2 vertical displacement', 'M2 angular displacement', 
                       'Source horizontal displacement', 'Source vertical displacement']

# for i in range(len(component_variation)):
#     calc_inputs = [0, 0, 0, 0, 0, 0, 0, 0, 0, False]
#     mean_distances = []

#     if 'angular' in component_variation[i]:
#         disp = np.concatenate((-np.flip(angular_displacement), angular_displacement))
#     else:
#         disp = np.concatenate((-np.flip(lin_displacement), lin_displacement))
    
#     for j in range(len(disp)):
#         calc_inputs[i] = disp[j]
#         mean_distances.append(setup_full_calc(*calc_inputs))
#     plt.scatter(disp, mean_distances)
#     plt.xlabel(f'{component_variation[i]} (deg))' if 'angular' in component_variation[i] else f'{component_variation[i]} (mm))')
#     plt.ylabel('Alignment FOM')
#     plt.xscale('symlog')
#     plt.savefig(f"C:/Users/etrem/Documents/Maitrise/Code/Results/OAP_sims/{component_variation[i]}", dpi=300, bbox_inches='tight')
#     plt.close()

calc_inputs = [0, 0, 0, 0, 0, 0, 0, 1e-2, 1e-2, True]
print(setup_full_calc(*calc_inputs))

# The precision of the code seems to be limitted by the intersection parameter in the find_beam_mirror_int fucntion. The lower the 
# intersection parameter, the more precise is the convergence of the beams under convergent conditions (focal, mirror face to face, etc.).