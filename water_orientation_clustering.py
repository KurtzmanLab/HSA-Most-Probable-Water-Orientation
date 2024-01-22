#!/usr/bin/python

# ===============================================================================
#
#          FILE: water_orientation_clustering.py
#
#         USAGE: ./water_orientation_clustering.py 'name_of_cluster.pdb'
#
#   DESCRIPTION: This code finds all probable orientations of water in each HSA cluster
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Mossa Ghattas, mossa.ghattas@gmail.com
#  ORGANIZATION: Lehman College, Chemistry Department
#       CREATED: 03/27/2020 11:10:42
#      REVISION:  ---
# ====================================================================

import numpy as np
import sys
import copy


####################################################################################
def perform_rotation_on_atom(q_0, q_1, q_2, q_3, atom_old_x_coords, atom_old_y_coords, atom_old_z_coords):

    atom_new_x_coords = np.multiply(np.subtract(np.add(np.square(q_0), np.square(q_1)), np.add(np.square(q_2), np.square(q_3))), atom_old_x_coords)
    atom_new_x_coords += np.multiply(2, np.multiply(np.add(np.multiply(q_1, q_2), np.multiply(q_0, q_3)), atom_old_y_coords))
    atom_new_x_coords += np.multiply(2, np.multiply(np.subtract(np.multiply(q_1, q_3), np.multiply(q_0, q_2)), atom_old_z_coords))

    atom_new_y_coords = np.multiply(np.multiply(2, np.subtract(np.multiply(q_1, q_2), np.multiply(q_0, q_3))), atom_old_x_coords)
    atom_new_y_coords += np.multiply(np.add(np.subtract(np.multiply(q_0, q_0), np.multiply(q_1, q_1)), np.subtract(np.multiply(q_2, q_2), np.multiply(q_3, q_3))), atom_old_y_coords)
    atom_new_y_coords += np.multiply(2, np.multiply(np.add(np.multiply(q_2, q_3), np.multiply(q_0, q_1)), atom_old_z_coords))

    atom_new_z_coords = np.multiply(np.multiply(2, np.add(np.multiply(q_1, q_3), np.multiply(q_0, q_2))), atom_old_x_coords)
    atom_new_z_coords += np.multiply(2, np.multiply(np.subtract(np.multiply(q_2, q_3), np.multiply(q_0, q_1)), atom_old_y_coords))
    atom_new_z_coords += np.multiply(np.add(np.subtract(np.subtract(np.multiply(q_0, q_0), np.multiply(q_1, q_1)), np.multiply(q_2, q_2)), np.multiply(q_3, q_3)), atom_old_z_coords)

    return [atom_new_x_coords, atom_new_y_coords, atom_new_z_coords]

####################################################################################


try:
    inputfile = sys.argv[1]

except IndexError:
    print("\nUSAGE: {} inputfile\n".format(sys.argv[0]))
    exit(2)

with open(inputfile, "r") as file:  # opens the file and reads each line as an item in a list and closes it.
    np.seterr(divide='ignore', invalid='ignore')  # https://stackoverflow.com/questions/14861891/runtimewarning-invalid-value-encountered-in-divide
    in_file = file.readlines()  # converts what is read to a list
    end = len(in_file)
    oxygen = in_file[1:end:3]  # extracting oxygen from in_file
    h1 = in_file[2:end:3]  # extracting h1 from in_file
    h2 = in_file[3:end:3]  # extracting h2 from in_file

####################################################################################

####################################################################################
#  Pulling the O atoms coordinates
oxygen_x = []
oxygen_y = []
oxygen_z = []
for i in oxygen:
    oxygen_x.append(float(i[32:38]))
    oxygen_y.append(float(i[40:46]))
    oxygen_z.append(float(i[48:54]))
# Turn the list of coordinates into a numpy array
np_o_x = np.array(oxygen_x)
np_o_y = np.array(oxygen_y)
np_o_z = np.array(oxygen_z)

####################################################################################

####################################################################################
# Pulling the H1 atoms coordinates
h1_x = []
h1_y = []
h1_z = []
for i in h1:
    h1_x.append(float(i[32:38]))
    h1_y.append(float(i[40:46]))
    h1_z.append(float(i[48:54]))
# Turn the list of coordinates into a numpy array
np_h1_x = np.array(h1_x)
np_h1_y = np.array(h1_y)
np_h1_z = np.array(h1_z)

####################################################################################

####################################################################################
# Pulling the H2 atoms coordinates
h2_x = []
h2_y = []
h2_z = []
for i in h2:
    h2_x.append(float(i[32:38]))
    h2_y.append(float(i[40:46]))
    h2_z.append(float(i[48:54]))
# Turn the list of coordinates into a numpy array
np_h2_x = np.array(h2_x)
np_h2_y = np.array(h2_y)
np_h2_z = np.array(h2_z)

####################################################################################

####################################################################################
# Translating all the H2O molecules such that all the O atoms are at the (0,0,0) position
# without altering the 3D orientation of the molecule
new_oxygen_x = np.subtract(np_o_x, np_o_x)
new_oxygen_y = np.subtract(np_o_y, np_o_y)
new_oxygen_z = np.subtract(np_o_z, np_o_z)

new_h1_x = np.subtract(np_h1_x, np_o_x)
new_h1_y = np.subtract(np_h1_y, np_o_y)
new_h1_z = np.subtract(np_h1_z, np_o_z)

new_h2_x = np.subtract(np_h2_x, np_o_x)
new_h2_y = np.subtract(np_h2_y, np_o_y)
new_h2_z = np.subtract(np_h2_z, np_o_z)

###############################################################################

###############################################################################
# Getting the length of h1 vector for each water molecule.
np.set_printoptions(suppress=True,
                    formatter={'float_kind': '{:0.8f}'.format})
h1_x_squares = np.square(new_h1_x)
h1_y_squares = np.square(new_h1_y)
h1_z_squares = np.square(new_h1_z)
h1_xy_squares_added = np.add(h1_x_squares, h1_y_squares)
h1_all_squares_added = np.add(h1_xy_squares_added, h1_z_squares)
h1_length = np.sqrt(h1_all_squares_added)

#  Getting the length of h2 vector for each water molecule.
np.set_printoptions(suppress=True,
                    formatter={'float_kind': '{:0.8f}'.format})
h2_x_squares = np.square(new_h2_x)
h2_y_squares = np.square(new_h2_y)
h2_z_squares = np.square(new_h2_z)
h2_xy_squares_added = np.add(h2_x_squares, h2_y_squares)
h2_all_squares_added = np.add(h2_xy_squares_added, h2_z_squares)
h2_length = np.sqrt(h2_all_squares_added)

###############################################################################

###############################################################################
# Normalizing each axis of each hydrogen atom
# We do this by dividing each atom axes with its own atom magnitude from oxygen atom
norm_h1_x = np.divide(new_h1_x, h1_length)
norm_h1_y = np.divide(new_h1_y, h1_length)
norm_h1_z = np.divide(new_h1_z, h1_length)

norm_h2_x = np.divide(new_h2_x, h2_length)
norm_h2_y = np.divide(new_h2_y, h2_length)
norm_h2_z = np.divide(new_h2_z, h2_length)

###############################################################################

###############################################################################
# Using this temporary which we will alter for all operations of rotations
temp_norm_h1_x = copy.deepcopy(norm_h1_x)
temp_norm_h1_y = copy.deepcopy(norm_h1_y)
temp_norm_h1_z = copy.deepcopy(norm_h1_z)
temp_norm_h2_x = copy.deepcopy(norm_h2_x)
temp_norm_h2_y = copy.deepcopy(norm_h2_y)
temp_norm_h2_z = copy.deepcopy(norm_h2_z)

###############################################################################

###############################################################################
# Defining the three axes of reference
x_ref = [1, 0, 0]
y_ref = [0, 1, 0]
z_ref = [0, 0, 1]

###############################################################################
# Rotating each water molecule such that: H1 is laid on the (1,0,0) AND H2 is laid planar to the xy plane.
H1final_q_0_all = []
H1final_q_1_all = []
H1final_q_2_all = []
H1final_q_3_all = []
H2final_q_0_all = []
H2final_q_1_all = []
H2final_q_2_all = []
H2final_q_3_all = []

H1updated_theta_in_degrees_all = []
H2updated_theta_in_degrees_all = []
H1updated_theta3p_in_degrees_all = []
H2updated_theta3p_in_degrees_all = []

for ndx, i in enumerate(norm_h1_x):
    #######################################################################################################################################################
    # This whole block is calculating first and second quaternions and their product quaternion to
    # ROTATE H1 on top of (1,0,0) and ROTATE H2 on the xy plane

    # Cross product of the two vectors will result in the axis of rotation - The first rotational axis (ar)
    # if O-H1 vector is parallel to x_ref vector, cross product is 0 and therefore we wont have an axis of rotation.
    # In the case where the O-H1 vector is lying on the (1,0,0) or the (-1,0,0) coordinates, we need to use y-ref (0,1,0) to just define cross product (the axis of rotation).
    # otherwise, we are good to continue using x_ref (1,0,0) for the cross product with the O-H1 vector
    if abs(temp_norm_h1_x[ndx]) == 1 and temp_norm_h1_y[ndx] == 0 and temp_norm_h1_z[ndx] == 0:
        H1ar_x = np.subtract(np.multiply(temp_norm_h1_y[ndx], y_ref[2]), np.multiply(temp_norm_h1_z[ndx], y_ref[1]))
        H1ar_y = np.subtract(np.multiply(temp_norm_h1_z[ndx], y_ref[0]), np.multiply(temp_norm_h1_x[ndx], y_ref[2]))
        H1ar_z = np.subtract(np.multiply(temp_norm_h1_x[ndx], y_ref[1]), np.multiply(temp_norm_h1_y[ndx], y_ref[0]))
    else:
        H1ar_x = np.subtract(np.multiply(temp_norm_h1_y[ndx], x_ref[2]), np.multiply(temp_norm_h1_z[ndx], x_ref[1]))
        H1ar_y = np.subtract(np.multiply(temp_norm_h1_z[ndx], x_ref[0]), np.multiply(temp_norm_h1_x[ndx], x_ref[2]))
        H1ar_z = np.subtract(np.multiply(temp_norm_h1_x[ndx], x_ref[1]), np.multiply(temp_norm_h1_y[ndx], x_ref[0]))

    # Normalizing the first rotational axis (ar).
    H1ar_x_squares = np.square(H1ar_x)
    H1ar_y_squares = np.square(H1ar_y)
    H1ar_z_squares = np.square(H1ar_z)
    H1ar_xy_squares_added = np.add(H1ar_x_squares, H1ar_y_squares)
    H1ar_all_squares_added = np.add(H1ar_xy_squares_added, H1ar_z_squares)
    H1ar_length = np.sqrt(H1ar_all_squares_added)

    H1norm_ar_x = np.nan_to_num(np.divide(H1ar_x, H1ar_length))
    H1norm_ar_y = np.nan_to_num(np.divide(H1ar_y, H1ar_length))
    H1norm_ar_z = np.nan_to_num(np.divide(H1ar_z, H1ar_length))

    # Calculating theta from the dot product of x_ref vector and h1 vector
    # https://math.stackexchange.com/questions/654315/how-to-convert-a-dot-product-of-two-vectors-to-the-angle-between-the-vectors
    H1dotproduct_h1_xref_xcomponent = np.multiply(temp_norm_h1_x[ndx], x_ref[0])
    H1dotproduct_h1_xref_ycomponent = np.multiply(temp_norm_h1_y[ndx], x_ref[1])  # I know it's zero and that's because x_ref is [1.0.0]
    H1dotproduct_h1_xref_zcomponent = np.multiply(temp_norm_h1_z[ndx], x_ref[2])  # I know it's zero and that's because x_ref is [1.0.0]
    H1dotproduct_h1_xref_xycomponents = np.add(H1dotproduct_h1_xref_xcomponent, H1dotproduct_h1_xref_ycomponent)
    H1dotproduct_h1_xref_allcomponents = np.add(H1dotproduct_h1_xref_xycomponents, H1dotproduct_h1_xref_zcomponent)
    H1theta_in_radians = np.arccos(H1dotproduct_h1_xref_allcomponents)
    H1theta_in_degrees = np.degrees(H1theta_in_radians)

    # Test to check for theta sign.
    H1dotproduct_ar_normh1_xcomponent = np.multiply(H1ar_x, temp_norm_h1_x[ndx])
    H1dotproduct_ar_normh1_ycomponent = np.multiply(H1ar_y, temp_norm_h1_y[ndx])
    H1dotproduct_ar_normh1_zcomponent = np.multiply(H1ar_z, temp_norm_h1_z[ndx])
    H1dotproduct_ar_normh1_xycomponents = np.add(H1dotproduct_ar_normh1_xcomponent, H1dotproduct_ar_normh1_ycomponent)
    H1dotproduct_ar_normh1_allcomponents = np.add(H1dotproduct_ar_normh1_xycomponents, H1dotproduct_ar_normh1_zcomponent)

    H1updated_theta_in_radians = np.zeros_like(H1theta_in_radians)
    if H1dotproduct_ar_normh1_allcomponents > 0:
        H1updated_theta_in_radians = H1theta_in_radians / 2.0
    else:
        H1updated_theta_in_radians = H1theta_in_radians / -2.0
    H1updated_theta_in_degrees = np.degrees(H1updated_theta_in_radians)
    H1updated_theta_in_degrees_all.append(H1updated_theta_in_degrees.tolist())

    # Defining first quaternion which will rotate the molecule such that H1 is on the (1,0,0) coordinates - based on the first axis of rotation (ar)
    H1q1_0 = np.cos(H1updated_theta_in_radians)
    H1q1_1 = np.multiply(H1norm_ar_x, np.sin(H1updated_theta_in_radians))
    H1q1_2 = np.multiply(H1norm_ar_y, np.sin(H1updated_theta_in_radians))
    H1q1_3 = np.multiply(H1norm_ar_z, np.sin(H1updated_theta_in_radians))

    e_0 = copy.deepcopy(H1q1_0)
    e_1 = copy.deepcopy(H1q1_1)
    e_2 = copy.deepcopy(H1q1_2)
    e_3 = copy.deepcopy(H1q1_3)

    # Performing that first quaternion on H1 and pull the new coordinates of H1
    new_h1_coords = perform_rotation_on_atom(e_0, e_1, e_2, e_3, temp_norm_h1_x[ndx], temp_norm_h1_y[ndx], temp_norm_h1_z[ndx])

    # H1 coordinates after the rotation operation - all the H1 coordinates should now be 1,0,0
    h1_moved_x = copy.deepcopy(new_h1_coords[0])
    h1_moved_y = copy.deepcopy(new_h1_coords[1])
    h1_moved_z = copy.deepcopy(new_h1_coords[2])

    # Performing that first quaternion on H2 and pull the new coordinates of H2
    new_h2_coords = perform_rotation_on_atom(e_0, e_1, e_2, e_3, temp_norm_h2_x[ndx], temp_norm_h2_y[ndx], temp_norm_h2_z[ndx])

    # H2 coordinates after the rotation operation - the H2 arent yet on the xy plane - the next pblock of code will take care of performing that
    h2_after_h1_moved_x = copy.deepcopy(new_h2_coords[0])
    h2_after_h1_moved_y = copy.deepcopy(new_h2_coords[1])
    h2_after_h1_moved_z = copy.deepcopy(new_h2_coords[2])

    # Cross product of the two vectors will result in the axis of rotation - The second rotational axis
    H1z_mol_vect_x = np.subtract(np.multiply(h1_moved_y, h2_after_h1_moved_z), np.multiply(h1_moved_z, h2_after_h1_moved_y))
    H1z_mol_vect_y = np.subtract(np.multiply(h1_moved_z, h2_after_h1_moved_x), np.multiply(h1_moved_x, h2_after_h1_moved_z))
    H1z_mol_vect_z = np.subtract(np.multiply(h1_moved_x, h2_after_h1_moved_y), np.multiply(h1_moved_y, h2_after_h1_moved_x))

    # Normalizing the second rotational axis.
    H1z_mol_vect_x_squares = np.square(H1z_mol_vect_x)
    H1z_mol_vect_y_squares = np.square(H1z_mol_vect_y)
    H1z_mol_vect_z_squares = np.square(H1z_mol_vect_z)
    H1z_mol_vect_xy_squares_added = np.add(H1z_mol_vect_x_squares, H1z_mol_vect_y_squares)
    H1z_mol_vect_all_squares_added = np.add(H1z_mol_vect_xy_squares_added, H1z_mol_vect_z_squares)
    H1z_mol_vect_length = np.sqrt(H1z_mol_vect_all_squares_added)

    H1norm_z_mol_vect_x = np.divide(H1z_mol_vect_x, H1z_mol_vect_length)
    H1norm_z_mol_vect_y = np.divide(H1z_mol_vect_y, H1z_mol_vect_length)
    H1norm_z_mol_vect_z = np.divide(H1z_mol_vect_z, H1z_mol_vect_length)

    # Calculating theta from the dot product of H1norm_z_mol_vect vector and z_ref vector
    # https://math.stackexchange.com/questions/654315/how-to-convert-a-dot-product-of-two-vectors-to-the-angle-between-the-vectors
    H1dotproduct_norm_z_mol_vect_zref_xcomponent = np.multiply(H1norm_z_mol_vect_x, z_ref[0])  # I know it's zero and that's because z_ref is [0,0,1]
    H1dotproduct_norm_z_mol_vect_zref_ycomponent = np.multiply(H1norm_z_mol_vect_y, z_ref[1])  # I know it's zero and that's because z_ref is [0,0,1]
    H1dotproduct_norm_z_mol_vect_zref_zcomponent = np.multiply(H1norm_z_mol_vect_z, z_ref[2])
    H1dotproduct_norm_z_mol_vect_zref_xycomponents = np.add(H1dotproduct_norm_z_mol_vect_zref_xcomponent, H1dotproduct_norm_z_mol_vect_zref_ycomponent)
    H1dotproduct_norm_z_mol_vect_zref_allcomponents = np.add(H1dotproduct_norm_z_mol_vect_zref_xycomponents, H1dotproduct_norm_z_mol_vect_zref_zcomponent)
    H1theta3p_in_radians = np.arccos(H1dotproduct_norm_z_mol_vect_zref_allcomponents)
    H1theta3p_in_degrees = np.degrees(H1theta3p_in_radians)

    # Test to check for theta sign.
    H1crossp_norm_z_mol_vect_z_ref_x = np.subtract(np.multiply(H1norm_z_mol_vect_y, z_ref[2]), np.multiply(H1norm_z_mol_vect_z, z_ref[1]))
    H1crossp_norm_z_mol_vect_z_ref_y = np.subtract(np.multiply(H1norm_z_mol_vect_z, z_ref[0]), np.multiply(H1norm_z_mol_vect_x, z_ref[2]))
    H1crossp_norm_z_mol_vect_z_ref_z = np.subtract(np.multiply(H1norm_z_mol_vect_x, z_ref[1]), np.multiply(H1norm_z_mol_vect_y, z_ref[0]))

    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent = np.multiply(H1crossp_norm_z_mol_vect_z_ref_x, h1_moved_x)
    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent = np.multiply(H1crossp_norm_z_mol_vect_z_ref_y, h1_moved_y)
    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent = np.multiply(H1crossp_norm_z_mol_vect_z_ref_z, h1_moved_z)
    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xycomponents = np.add(H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent, H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent)
    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents = np.add(H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xycomponents, H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent)

    H1updated_theta3p_in_radians = np.zeros_like(H1theta3p_in_radians)

    if H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents < 0.0:
        H1updated_theta3p_in_radians = H1theta3p_in_radians / 2.0
    else:
        H1updated_theta3p_in_radians = H1theta3p_in_radians / -2.0

    H1updated_theta3p_in_degrees = np.degrees(H1updated_theta3p_in_radians)
    H1updated_theta3p_in_degrees_all.append(H1updated_theta3p_in_degrees.tolist())

    # Defining Second quaternion which will rotate the molecule such that H1 is STILL on the (1,0,0) coordinates but now the H2 will be rotated such that it lies planar to the xy plane
    H1q2_0 = np.cos(H1updated_theta3p_in_radians)
    H1q2_1 = np.multiply(x_ref[0], np.sin(H1updated_theta3p_in_radians))
    H1q2_2 = np.multiply(x_ref[1], np.sin(H1updated_theta3p_in_radians))
    H1q2_3 = np.multiply(x_ref[2], np.sin(H1updated_theta3p_in_radians))

    #  Calculating the product quaternion.
    #  http://sneg.co.uk/content/quick-guide-quaternions
    #  http://sneg.co.uk/content/quick-guide-quaternions
    #  http://answers.google.com/answers/threadview/id/596035.html
    #  http://en.wikipedia.org/wiki/Quaternion#Multiplication_of_basis_elements
    #  http://www.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation#Multiplying_quaternions
    #  http://wiki.alioth.net/index.php/Quaternion
    #  http://mathworld.wolfram.com/Quaternion.html
    #  http://www.csse.uwa.edu.au/~pk/research/matlabfns/Rotations/quaternionproduct.m
    #  http://www.csse.uwa.edu.au/~pk/research/matlabfns/Rotations/quaternionproduct.m
    e_0 = np.subtract(np.subtract(np.subtract(np.multiply(H1q1_0, H1q2_0), np.multiply(H1q1_1, H1q2_1)), np.multiply(H1q1_2, H1q2_2)), np.multiply(H1q1_3, H1q2_3))
    e_1 = np.subtract(np.add(np.add(np.multiply(H1q1_0, H1q2_1), np.multiply(H1q1_1, H1q2_0)), np.multiply(H1q1_2, H1q2_3)), np.multiply(H1q1_3, H1q2_2))
    e_2 = np.add(np.add(np.subtract(np.multiply(H1q1_0, H1q2_2), np.multiply(H1q1_1, H1q2_3)), np.multiply(H1q1_2, H1q2_0)), np.multiply(H1q1_3, H1q2_1))
    e_3 = np.add(np.subtract(np.add(np.multiply(H1q1_0, H1q2_3), np.multiply(H1q1_1, H1q2_2)), np.multiply(H1q1_2, H1q2_1)), np.multiply(H1q1_3, H1q2_0))

    H1final_q_0 = copy.deepcopy(e_0)
    H1final_q_1 = copy.deepcopy(e_1)
    H1final_q_2 = copy.deepcopy(e_2)
    H1final_q_3 = copy.deepcopy(e_3)

    H1final_q_0_all.append(H1final_q_0)
    H1final_q_1_all.append(H1final_q_1)
    H1final_q_2_all.append(H1final_q_2)
    H1final_q_3_all.append(H1final_q_3)
    # This whole block was calculating first and second quaternions and their product quaternion that will rotate H1 on top of 1,0,0 and H2 on the xy plane
    #######################################################################################################################################################

    #######################################################################################################################################################
    # This whole block is calculating first and second quaternions and their product quaternion to
    # ROTATE H2 on top of (1,0,0) and ROTATE H1 on the xy plane

    # Cross product of the two vectors will result in the axis of rotation - The first rotational axis (ar)
    # if O-H2 vector is parallel to x_ref vector, cross product is 0 and therefore we wont have an axis of rotation.
    # In the case where the O-H2 vector is lying on the (1,0,0) or the (-1,0,0) coordinates, we need to use y-ref (0,1,0) to just define cross product (the axis of rotation).
    # otherwise, we are good to continue using x_ref (1,0,0) for the cross product with the O-H1 vector
    if abs(temp_norm_h2_x[ndx]) == 1 and temp_norm_h2_y[ndx] == 0 and temp_norm_h2_z[ndx] == 0:
        H2ar_x = np.subtract(np.multiply(temp_norm_h2_y[ndx], y_ref[2]), np.multiply(temp_norm_h2_z[ndx], y_ref[1]))
        H2ar_y = np.subtract(np.multiply(temp_norm_h2_z[ndx], y_ref[0]), np.multiply(temp_norm_h2_x[ndx], y_ref[2]))
        H2ar_z = np.subtract(np.multiply(temp_norm_h2_x[ndx], y_ref[1]), np.multiply(temp_norm_h2_y[ndx], y_ref[0]))
    else:
        H2ar_x = np.subtract(np.multiply(temp_norm_h2_y[ndx], x_ref[2]), np.multiply(temp_norm_h2_z[ndx], x_ref[1]))
        H2ar_y = np.subtract(np.multiply(temp_norm_h2_z[ndx], x_ref[0]), np.multiply(temp_norm_h2_x[ndx], x_ref[2]))
        H2ar_z = np.subtract(np.multiply(temp_norm_h2_x[ndx], x_ref[1]), np.multiply(temp_norm_h2_y[ndx], x_ref[0]))

    #  Normalizing the first rotational axis (ar).
    H2ar_x_squares = np.square(H2ar_x)
    H2ar_y_squares = np.square(H2ar_y)
    H2ar_z_squares = np.square(H2ar_z)
    H2ar_xy_squares_added = np.add(H2ar_x_squares, H2ar_y_squares)
    H2ar_all_squares_added = np.add(H2ar_xy_squares_added, H2ar_z_squares)
    H2ar_length = np.sqrt(H2ar_all_squares_added)

    H2norm_ar_x = np.nan_to_num(np.divide(H2ar_x, H2ar_length))
    H2norm_ar_y = np.nan_to_num(np.divide(H2ar_y, H2ar_length))
    H2norm_ar_z = np.nan_to_num(np.divide(H2ar_z, H2ar_length))

    #  Calculating theta from the dot product of x_ref vector and h1 vector
    #  https://math.stackexchange.com/questions/654315/how-to-convert-a-dot-product-of-two-vectors-to-the-angle-between-the-vectors
    H2dotproduct_h2_xref_xcomponent = np.multiply(temp_norm_h2_x[ndx], x_ref[0])
    H2dotproduct_h2_xref_ycomponent = np.multiply(temp_norm_h2_y[ndx], x_ref[1])  # I know it's zero and that's because x_ref is [1.0.0]
    H2dotproduct_h2_xref_zcomponent = np.multiply(temp_norm_h2_z[ndx], x_ref[2])  # I know it's zero and that's because x_ref is [1.0.0]
    H2dotproduct_h2_xref_xycomponents = np.add(H2dotproduct_h2_xref_xcomponent, H2dotproduct_h2_xref_ycomponent)
    H2dotproduct_h2_xref_allcomponents = np.add(H2dotproduct_h2_xref_xycomponents, H2dotproduct_h2_xref_zcomponent)
    H2theta_in_radians = np.arccos(H2dotproduct_h2_xref_allcomponents)
    H2theta_in_degrees = np.degrees(H2theta_in_radians)

    # Test to check for theta sign.
    H2dotproduct_ar_normh2_xcomponent = np.multiply(H2ar_x, temp_norm_h2_x[ndx])
    H2dotproduct_ar_normh2_ycomponent = np.multiply(H2ar_y, temp_norm_h2_y[ndx])
    H2dotproduct_ar_normh2_zcomponent = np.multiply(H2ar_z, temp_norm_h2_z[ndx])
    H2dotproduct_ar_normh2_xycomponents = np.add(H2dotproduct_ar_normh2_xcomponent, H2dotproduct_ar_normh2_ycomponent)
    H2dotproduct_ar_normh2_allcomponents = np.add(H2dotproduct_ar_normh2_xycomponents, H2dotproduct_ar_normh2_zcomponent)

    H2updated_theta_in_radians = np.zeros_like(H2theta_in_radians)

    if H2dotproduct_ar_normh2_allcomponents > 0:
        H2updated_theta_in_radians = H2theta_in_radians / 2.0
    else:
        H2updated_theta_in_radians = H2theta_in_radians / -2.0

    H2updated_theta_in_degrees = np.degrees(H2updated_theta_in_radians)
    H2updated_theta_in_degrees_all.append(H2updated_theta_in_degrees.tolist())

    # Defining first quaternion based on the first axis of rotation (ar).
    H2q1_0 = np.cos(H2updated_theta_in_radians)
    H2q1_1 = np.multiply(H2norm_ar_x, np.sin(H2updated_theta_in_radians))
    H2q1_2 = np.multiply(H2norm_ar_y, np.sin(H2updated_theta_in_radians))
    H2q1_3 = np.multiply(H2norm_ar_z, np.sin(H2updated_theta_in_radians))

    e_0 = copy.deepcopy(H2q1_0)
    e_1 = copy.deepcopy(H2q1_1)
    e_2 = copy.deepcopy(H2q1_2)
    e_3 = copy.deepcopy(H2q1_3)

    # Performing that first quaternion on H2 and pull the new coordinates of H2
    new_h2_coords = perform_rotation_on_atom(e_0, e_1, e_2, e_3, temp_norm_h2_x[ndx], temp_norm_h2_y[ndx], temp_norm_h2_z[ndx])

    # H2 coordinates after the rotation operation - all the H2 coordinates should now be 1,0,0
    h2_moved_x = copy.deepcopy(new_h2_coords[0])
    h2_moved_y = copy.deepcopy(new_h2_coords[1])
    h2_moved_z = copy.deepcopy(new_h2_coords[2])

    # Performing that first quaternion on H1 and pull the new coordinates of H1
    new_h1_coords = perform_rotation_on_atom(e_0, e_1, e_2, e_3, temp_norm_h1_x[ndx], temp_norm_h1_y[ndx], temp_norm_h1_z[ndx])

    # H1 coordinates after the rotation operation - the H1 aren't yet on the xy plane - the next block of code will take care of performing that
    h1_after_h2_moved_x = copy.deepcopy(new_h1_coords[0])
    h1_after_h2_moved_y = copy.deepcopy(new_h1_coords[1])
    h1_after_h2_moved_z = copy.deepcopy(new_h1_coords[2])

    # Cross product of the two vectors will result in the axis of rotation - The second rotational axis
    H2z_mol_vect_x = np.subtract(np.multiply(h2_moved_y, h1_after_h2_moved_z), np.multiply(h2_moved_z, h1_after_h2_moved_y))
    H2z_mol_vect_y = np.subtract(np.multiply(h2_moved_z, h1_after_h2_moved_x), np.multiply(h2_moved_x, h1_after_h2_moved_z))
    H2z_mol_vect_z = np.subtract(np.multiply(h2_moved_x, h1_after_h2_moved_y), np.multiply(h2_moved_y, h1_after_h2_moved_x))

    # Normalizing the second rotational axis.
    H2z_mol_vect_x_squares = np.square(H2z_mol_vect_x)
    H2z_mol_vect_y_squares = np.square(H2z_mol_vect_y)
    H2z_mol_vect_z_squares = np.square(H2z_mol_vect_z)
    H2z_mol_vect_xy_squares_added = np.add(H2z_mol_vect_x_squares, H2z_mol_vect_y_squares)
    H2z_mol_vect_all_squares_added = np.add(H2z_mol_vect_xy_squares_added, H2z_mol_vect_z_squares)
    H2z_mol_vect_length = np.sqrt(H2z_mol_vect_all_squares_added)

    H2norm_z_mol_vect_x = np.divide(H2z_mol_vect_x, H2z_mol_vect_length)
    H2norm_z_mol_vect_y = np.divide(H2z_mol_vect_y, H2z_mol_vect_length)
    H2norm_z_mol_vect_z = np.divide(H2z_mol_vect_z, H2z_mol_vect_length)

    # Calculating theta from the dot product of H1norm_z_mol_vect vector and z_ref vector
    # https://math.stackexchange.com/questions/654315/how-to-convert-a-dot-product-of-two-vectors-to-the-angle-between-the-vectors
    H2dotproduct_norm_z_mol_vect_zref_xcomponent = np.multiply(H2norm_z_mol_vect_x, z_ref[0])  # I know it's zero and that's because z_ref is [0,0,1]
    H2dotproduct_norm_z_mol_vect_zref_ycomponent = np.multiply(H2norm_z_mol_vect_y, z_ref[1])  # I know it's zero and that's because z_ref is [0,0,1]
    H2dotproduct_norm_z_mol_vect_zref_zcomponent = np.multiply(H2norm_z_mol_vect_z, z_ref[2])
    H2dotproduct_norm_z_mol_vect_zref_xycomponents = np.add(H2dotproduct_norm_z_mol_vect_zref_xcomponent, H2dotproduct_norm_z_mol_vect_zref_ycomponent)
    H2dotproduct_norm_z_mol_vect_zref_allcomponents = np.add(H2dotproduct_norm_z_mol_vect_zref_xycomponents, H2dotproduct_norm_z_mol_vect_zref_zcomponent)

    H2theta3p_in_radians = np.arccos(H2dotproduct_norm_z_mol_vect_zref_allcomponents)
    H2theta3p_in_degrees = np.degrees(H2theta3p_in_radians)

    # Test to check for theta sign.
    H2crossp_norm_z_mol_vect_z_ref_x = np.subtract(np.multiply(H2norm_z_mol_vect_y, z_ref[2]), np.multiply(H2norm_z_mol_vect_z, z_ref[1]))
    H2crossp_norm_z_mol_vect_z_ref_y = np.subtract(np.multiply(H2norm_z_mol_vect_z, z_ref[0]), np.multiply(H2norm_z_mol_vect_x, z_ref[2]))
    H2crossp_norm_z_mol_vect_z_ref_z = np.subtract(np.multiply(H2norm_z_mol_vect_x, z_ref[1]), np.multiply(H2norm_z_mol_vect_y, z_ref[0]))

    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent = np.multiply(H2crossp_norm_z_mol_vect_z_ref_x, h2_moved_x)
    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent = np.multiply(H2crossp_norm_z_mol_vect_z_ref_y, h2_moved_y)
    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent = np.multiply(H2crossp_norm_z_mol_vect_z_ref_z, h2_moved_z)
    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xycomponents = np.add(H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent, H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent)
    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents = np.add(H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xycomponents, H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent)

    H2updated_theta3p_in_radians = np.zeros_like(H2theta3p_in_radians)

    if H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents < 0.0:
        H2updated_theta3p_in_radians = H2theta3p_in_radians / 2.0
    else:
        H2updated_theta3p_in_radians = H2theta3p_in_radians / -2.0

    H2updated_theta3p_in_degrees = np.degrees(H2updated_theta3p_in_radians)
    H2updated_theta3p_in_degrees_all.append(H2updated_theta3p_in_degrees.tolist())
    # Defining Second quaternion which will rotate the molecule such that H2 is STILL on the (1,0,0) coordinates but now the H1 will be rotated such that it lies planar to the xy plane
    H2q2_0 = np.cos(H2updated_theta3p_in_radians)
    H2q2_1 = np.multiply(x_ref[0], np.sin(H2updated_theta3p_in_radians))
    H2q2_2 = np.multiply(x_ref[1], np.sin(H2updated_theta3p_in_radians))
    H2q2_3 = np.multiply(x_ref[2], np.sin(H2updated_theta3p_in_radians))

    #  Calculating the product quaternion.
    #  http://sneg.co.uk/content/quick-guide-quaternions
    #  http://sneg.co.uk/content/quick-guide-quaternions
    #  http://answers.google.com/answers/threadview/id/596035.html
    #  http://en.wikipedia.org/wiki/Quaternion#Multiplication_of_basis_elements
    #  http://www.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation#Multiplying_quaternions
    #  http://wiki.alioth.net/index.php/Quaternion
    #  http://mathworld.wolfram.com/Quaternion.html
    #  http://www.csse.uwa.edu.au/~pk/research/matlabfns/Rotations/quaternionproduct.m
    #  http://www.csse.uwa.edu.au/~pk/research/matlabfns/Rotations/quaternionproduct.m
    e_0 = np.subtract(np.subtract(np.subtract(np.multiply(H2q1_0, H2q2_0), np.multiply(H2q1_1, H2q2_1)), np.multiply(H2q1_2, H2q2_2)), np.multiply(H2q1_3, H2q2_3))
    e_1 = np.subtract(np.add(np.add(np.multiply(H2q1_0, H2q2_1), np.multiply(H2q1_1, H2q2_0)), np.multiply(H2q1_2, H2q2_3)), np.multiply(H2q1_3, H2q2_2))
    e_2 = np.add(np.add(np.subtract(np.multiply(H2q1_0, H2q2_2), np.multiply(H2q1_1, H2q2_3)), np.multiply(H2q1_2, H2q2_0)), np.multiply(H2q1_3, H2q2_1))
    e_3 = np.add(np.subtract(np.add(np.multiply(H2q1_0, H2q2_3), np.multiply(H2q1_1, H2q2_2)), np.multiply(H2q1_2, H2q2_1)), np.multiply(H2q1_3, H2q2_0))

    H2final_q_0 = copy.deepcopy(e_0)
    H2final_q_1 = copy.deepcopy(e_1)
    H2final_q_2 = copy.deepcopy(e_2)
    H2final_q_3 = copy.deepcopy(e_3)

    H2final_q_0_all.append(H2final_q_0)
    H2final_q_1_all.append(H2final_q_1)
    H2final_q_2_all.append(H2final_q_2)
    H2final_q_3_all.append(H2final_q_3)
    # This whole block was calculating first and second quaternions and their product quaternion that will rotate H2 on top of 1,0,0 and H1 on the xy plane
    #######################################################################################################################################################

##########################################################################################################
# Initially the list to create a list of list with size N
N = len(H1final_q_0_all)
cutoff = 0.15000
Neighbors = [None] * N
for i in range(N):
    Neighbors[i] = []

for i in range(0, N):
    # clear_output(wait=True)  # to show in console the progress
    for j in range(i+1, N):
        if j != i:
            q_0_diff_AA = H1final_q_0_all[i] - H1final_q_0_all[j]
            q_1_diff_AA = H1final_q_1_all[i] - H1final_q_1_all[j]
            q_2_diff_AA = H1final_q_2_all[i] - H1final_q_2_all[j]
            q_3_diff_AA = H1final_q_3_all[i] - H1final_q_3_all[j]
            q_0_add_AA = H1final_q_0_all[i] + H1final_q_0_all[j]
            q_1_add_AA = H1final_q_1_all[i] + H1final_q_1_all[j]
            q_2_add_AA = H1final_q_2_all[i] + H1final_q_2_all[j]
            q_3_add_AA = H1final_q_3_all[i] + H1final_q_3_all[j]

            q_0_diff_AB = H1final_q_0_all[i] - H2final_q_0_all[j]
            q_1_diff_AB = H1final_q_1_all[i] - H2final_q_1_all[j]
            q_2_diff_AB = H1final_q_2_all[i] - H2final_q_2_all[j]
            q_3_diff_AB = H1final_q_3_all[i] - H2final_q_3_all[j]
            q_0_add_AB = H1final_q_0_all[i] + H2final_q_0_all[j]
            q_1_add_AB = H1final_q_1_all[i] + H2final_q_1_all[j]
            q_2_add_AB = H1final_q_2_all[i] + H2final_q_2_all[j]
            q_3_add_AB = H1final_q_3_all[i] + H2final_q_3_all[j]

            q_0_diff_BA = H2final_q_0_all[i] - H1final_q_0_all[j]
            q_1_diff_BA = H2final_q_1_all[i] - H1final_q_1_all[j]
            q_2_diff_BA = H2final_q_2_all[i] - H1final_q_2_all[j]
            q_3_diff_BA = H2final_q_3_all[i] - H1final_q_3_all[j]
            q_0_add_BA = H2final_q_0_all[i] + H1final_q_0_all[j]
            q_1_add_BA = H2final_q_1_all[i] + H1final_q_1_all[j]
            q_2_add_BA = H2final_q_2_all[i] + H1final_q_2_all[j]
            q_3_add_BA = H2final_q_3_all[i] + H1final_q_3_all[j]

            q_0_diff_BB = H2final_q_0_all[i] - H2final_q_0_all[j]
            q_1_diff_BB = H2final_q_1_all[i] - H2final_q_1_all[j]
            q_2_diff_BB = H2final_q_2_all[i] - H2final_q_2_all[j]
            q_3_diff_BB = H2final_q_3_all[i] - H2final_q_3_all[j]
            q_0_add_BB = H2final_q_0_all[i] + H2final_q_0_all[j]
            q_1_add_BB = H2final_q_1_all[i] + H2final_q_1_all[j]
            q_2_add_BB = H2final_q_2_all[i] + H2final_q_2_all[j]
            q_3_add_BB = H2final_q_3_all[i] + H2final_q_3_all[j]

            Euclidean_distance_by_q = min((2 * min(((q_0_diff_AA ** 2) + (q_1_diff_AA ** 2) + (q_2_diff_AA ** 2) + (q_3_diff_AA ** 2)), ((q_0_add_AA ** 2) + (q_1_add_AA ** 2) + (q_2_add_AA ** 2) + (q_3_add_AA ** 2)))), (2 * min(((q_0_diff_AB ** 2) + (q_1_diff_AB ** 2) + (q_2_diff_AB ** 2) + (q_3_diff_AB ** 2)), ((q_0_add_AB ** 2) + (q_1_add_AB ** 2) + (q_2_add_AB ** 2) + (q_3_add_AB ** 2)))) , (2 * min(((q_0_diff_BB ** 2) + (q_1_diff_BB ** 2) + (q_2_diff_BB ** 2) + (q_3_diff_BB ** 2)), ((q_0_add_BB ** 2) + (q_1_add_BB ** 2) + (q_2_add_BB ** 2) + (q_3_add_BB ** 2)))) , (2 * min(((q_0_diff_BA ** 2) + (q_1_diff_BA ** 2) + (q_2_diff_BA ** 2) + (q_3_diff_BA ** 2)), ((q_0_add_BA ** 2) + (q_1_add_BA ** 2) + (q_2_add_BA ** 2) + (q_3_add_BA ** 2)))))
            if Euclidean_distance_by_q < cutoff:
                Neighbors[i].append(j)
                Neighbors[j].append(i)


sorted_tuples = [None] * N
for i in range(N):
    sorted_tuples[i] = (len(Neighbors[i]), i)
##########################################################################################################


sorted_tuples.sort(reverse=True)

available = [1] * len(H1final_q_0_all)

most_neighbored_molecules = []  # most_neighbored_molecules is a list that I extracted from sorted_tuples to be able to handle better. most_neighbored_molecules contains H2O ID that has the most amount of neighbors.
for x in sorted_tuples:
    most_neighbored_molecules.append(x[1])

for i in most_neighbored_molecules:
    if available[i] == 1:
        j = copy.deepcopy(Neighbors[i])
        available[i] = 0
        for k in j:
            if available[k] == 1:
                available[k] = 0
                Neighbors[k] = []
            else:
                Neighbors[i].remove(k)
        if Neighbors[i] != []:
            Neighbors[i].append(i)
    else:
        Neighbors[i] = []

Neighbors_count = 0
for i in Neighbors:
    for j in i:
        Neighbors_count += 1

clusters = []
for ndx, i in enumerate(Neighbors):
    if len(i) != 0:
        clusters.append((ndx, len(i), round((len(i) / Neighbors_count) * 100, 3)))

sorted_cluster = copy.deepcopy(sorted(clusters, key=lambda x: x[2], reverse=True))

##########################################################################

#  This block is to test if the orientational clusters aren't insignificant.

sum_percenatge_of_clusters = 0
for ndx, i in enumerate(sorted_cluster):
    if i[2] >= 10.000:
        sum_percenatge_of_clusters += i[2]

final_orientational_clusters = []
if sum_percenatge_of_clusters >= 40:
    for ndx1, j in enumerate(sorted_cluster):
        if j[2] >= 10.000:
            final_orientational_clusters.append(j[0])
else:
    for ndx2, k in enumerate(sorted_cluster):
        if ndx2 <= 4:
            final_orientational_clusters.append(k[0])

###############################################################################

outfile = open("Water_orientational_clusters_of_{}".format(inputfile[0:18]), "w")
outfile.write("REMARKS\n")
outfile.write("This pdb file has the most probable orientations in a water cluster.\n")
outfile.write("We consider an orientational cluster if it has at least 10% of water population as neighbors.\n")
atmndx = 1
for ndx, i in enumerate(sorted_cluster):
    if i[0] in final_orientational_clusters:
        outfile.write("{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}  {}%\n".format("ATOM", atmndx, "OW", "SOL", "C", 5688, oxygen_x[i[0]], oxygen_y[i[0]], oxygen_z[i[0]], i[2]))
        atmndx += 1
        outfile.write("{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format("ATOM", atmndx, "H1", "SOL", "X", 5688, h1_x[i[0]], h1_y[i[0]], h1_z[i[0]]))
        atmndx += 1
        outfile.write("{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format("ATOM", atmndx, "H2", "SOL", "X", 5688, h2_x[i[0]], h2_y[i[0]], h2_z[i[0]]))
        atmndx += 1
outfile.close()

###############################################################################
