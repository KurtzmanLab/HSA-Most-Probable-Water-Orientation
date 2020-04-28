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

import math
import numpy as np
import sys
import copy
import datetime
import os

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

###############################################################################

#  Setting new coordinates for O for each water molecule.

oxygen_x = []
for i in oxygen:
    oxygen_x.append(float(i[32:38]))
np_o_x = np.array(oxygen_x)
#  print("oxygen_x:", np_o_x)
new_oxygen_x = np.subtract(np_o_x, np_o_x)
#  print("new_oxygen_x:", new_oxygen_x)

oxygen_y = []
for j in oxygen:
    oxygen_y.append(float(j[40:46]))
np_o_y = np.array(oxygen_y)
#  print("oxygen_y:", np_o_y)
new_oxygen_y = np.subtract(np_o_y, np_o_y)
#  print("new_oxygen_y:", new_oxygen_y)

oxygen_z = []
for k in oxygen:
    oxygen_z.append(float(k[48:54]))
np_o_z = np.array(oxygen_z)
#  print("oxygen_z:", np_o_z)
new_oxygen_z = np.subtract(np_o_z, np_o_z)
#  print("new_oxygen_z:", new_oxygen_z)
############

# Setting new coordinates for h1 with respect to O atom for each water molecule
h1_x = []
for i in h1:
    h1_x.append(float(i[32:38]))
np_h1_x = np.array(h1_x)
#  print("h1_x", np_h1_x)
new_h1_x = np.subtract(np_h1_x, np_o_x)
#  print("new_h1_x:", new_h1_x)

h1_y = []
for j in h1:
    h1_y.append(float(j[40:46]))
np_h1_y = np.array(h1_y)
#  print("h1_y", np_h1_y)
new_h1_y = np.subtract(np_h1_y, np_o_y)
#  print("new_h1_y:", new_h1_y)

h1_z = []
for k in h1:
    h1_z.append(float(k[48:54]))
np_h1_z = np.array(h1_z)
#  print("h1_z", np_h1_z)
new_h1_z = np.subtract(np_h1_z, np_o_z)
#  print("new_h1_z:", new_h1_z)
############

# Setting new coordinates for h2 with respect to O atom for each water molecule.
h2_x = []
for i in h2:
    h2_x.append(float(i[32:38]))
np_h2_x = np.array(h2_x)
#  print("h2_x:", np_h2_x)
new_h2_x = np.subtract(np_h2_x, np_o_x)
#  print("new_h2_x:", new_h2_x)

h2_y = []
for j in h2:
    h2_y.append(float(j[40:46]))
np_h2_y = np.array(h2_y)
#  print("h2_y:", np_h2_y)
new_h2_y = np.subtract(np_h2_y, np_o_y)
#  print("new_h2_y:", new_h2_y)

h2_z = []
for k in h2:
    h2_z.append(float(k[48:54]))
np_h2_z = np.array(h2_z)
#  print("h2_z:", np_h2_z)
new_h2_z = np.subtract(np_h2_z, np_o_z)
#  print("new_h2_z:", new_h2_z)

###############################################################################

#  Getting the length of h1 vector for each water molecule.

np.set_printoptions(suppress=True,
                    formatter={'float_kind': '{:0.8f}'.format})
h1_x_squares = np.square(new_h1_x)
#  print(h1_x_squares)
h1_y_squares = np.square(new_h1_y)
#  print(h1_y_squares)
h1_z_squares = np.square(new_h1_z)
#  print(h1_z_squares)
h1_xy_squares_added = np.add(h1_x_squares, h1_y_squares)
h1_all_squares_added = np.add(h1_xy_squares_added, h1_z_squares)
#  print(h1_all_squares_added)
h1_length = np.sqrt(h1_all_squares_added)
#  print("h1_length:", h1_length)


#  Getting the length of h2 vector for each water molecule.
np.set_printoptions(suppress=True,
                    formatter={'float_kind': '{:0.8f}'.format})
h2_x_squares = np.square(new_h2_x)
#  print(h2_x_squares)
h2_y_squares = np.square(new_h2_y)
#  print(h2_y_squares)
h2_z_squares = np.square(new_h2_z)
#  print(h2_z_squares)
h2_xy_squares_added = np.add(h2_x_squares, h2_y_squares)
h2_all_squares_added = np.add(h2_xy_squares_added, h2_z_squares)
#  print(h2_all_squares_added)
h2_length = np.sqrt(h2_all_squares_added)
#  print("h2_length:", h2_length)

###############################################################################

#  Normalizing each axis of each hydrogen atom. we do this by dividing each atom axes with its own atom magnitude from oxygen atom

norm_h1_x = np.divide(new_h1_x, h1_length)
#  print("norm_h1_x:", norm_h1_x)

norm_h1_y = np.divide(new_h1_y, h1_length)
#  print("norm_h1_y:", norm_h1_y)

norm_h1_z = np.divide(new_h1_z, h1_length)
#  print("norm_h1_z:", norm_h1_z)

norm_h2_x = np.divide(new_h2_x, h2_length)
#  print("norm_h2_x:", norm_h2_x)

norm_h2_y = np.divide(new_h2_y, h2_length)
#  print("norm_h2_y:", norm_h2_y)

norm_h2_z = np.divide(new_h2_z, h2_length)
#  print("norm_h2_z:", norm_h2_z)


x_ref = [1, 0, 0]
y_ref = [0, 1, 0]
z_ref = [0, 0, 1]

########################################################################################################################################

temp_norm_h1_x = copy.deepcopy(norm_h1_x)
temp_norm_h1_y = copy.deepcopy(norm_h1_y)
temp_norm_h1_z = copy.deepcopy(norm_h1_z)
temp_norm_h2_x = copy.deepcopy(norm_h2_x)
temp_norm_h2_y = copy.deepcopy(norm_h2_y)
temp_norm_h2_z = copy.deepcopy(norm_h2_z)

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

    H1ar_x = np.subtract(np.multiply(temp_norm_h1_y[ndx], x_ref[2]), np.multiply(temp_norm_h1_z[ndx], x_ref[1]))
    #  print("ar_x:", ar_x)
    H1ar_y = np.subtract(np.multiply(temp_norm_h1_z[ndx], x_ref[0]), np.multiply(temp_norm_h1_x[ndx], x_ref[2]))
    #  print("ar_y:", ar_y)
    H1ar_z = np.subtract(np.multiply(temp_norm_h1_x[ndx], x_ref[1]), np.multiply(temp_norm_h1_y[ndx], x_ref[0]))
    #  print("ar_z:", ar_z)

    #  Normalizing the first rotational axis (ar).

    H1ar_x_squares = np.square(H1ar_x)
    #  print(ar_x_squares)
    H1ar_y_squares = np.square(H1ar_y)
    #  print(ar_y_squares)
    H1ar_z_squares = np.square(H1ar_z)
    #  print(ar_z_squares)
    H1ar_xy_squares_added = np.add(H1ar_x_squares, H1ar_y_squares)
    H1ar_all_squares_added = np.add(H1ar_xy_squares_added, H1ar_z_squares)
    #  print(ar_all_squares_added)
    H1ar_length = np.sqrt(H1ar_all_squares_added)
    #  print("ar_length:", ar_length)

    H1norm_ar_x = np.nan_to_num(np.divide(H1ar_x, H1ar_length))
    #  print("norm_ar_x:", norm_ar_x)
    H1norm_ar_y = np.nan_to_num(np.divide(H1ar_y, H1ar_length))
    #  print("norm_ar_y:", norm_ar_y)
    H1norm_ar_z = np.nan_to_num(np.divide(H1ar_z, H1ar_length))
    #  print("norm_ar_z:", norm_ar_z)

    ###############################################################################

    #  Calculating theta from the dot product of x_ref vector and h1 vector
    #  https://math.stackexchange.com/questions/654315/how-to-convert-a-dot-product-of-two-vectors-to-the-angle-between-the-vectors

    H1dotproduct_h1_xref_xcomponent = np.multiply(temp_norm_h1_x[ndx], x_ref[0])
    H1dotproduct_h1_xref_ycomponent = np.multiply(temp_norm_h1_y[ndx], x_ref[1])  # I know it's zero and that's because x_ref is [1.0.0]
    H1dotproduct_h1_xref_zcomponent = np.multiply(temp_norm_h1_z[ndx], x_ref[2])  # I know it's zero and that's because x_ref is [1.0.0]
    H1dotproduct_h1_xref_xycomponents = np.add(H1dotproduct_h1_xref_xcomponent, H1dotproduct_h1_xref_ycomponent)
    H1dotproduct_h1_xref_allcomponents = np.add(H1dotproduct_h1_xref_xycomponents, H1dotproduct_h1_xref_zcomponent)
    #  print(dotproduct_h1_xref_allcomponents)
    H1theta_in_radians = np.arccos(H1dotproduct_h1_xref_allcomponents)
    #  print("theta_in_radians:", theta_in_radians)
    H1theta_in_degrees = np.degrees(H1theta_in_radians)
    #  print(theta_in_degrees)
    #  print(max(theta_in_degrees))
    #  print(min(theta_in_degrees))

    ###############################################################################

    #  Test to check for theta sign.

    H1dotproduct_ar_normh1_xcomponent = np.multiply(H1ar_x, temp_norm_h1_x[ndx])
    #  print(dotproduct_ar_normh1_xcomponent)
    H1dotproduct_ar_normh1_ycomponent = np.multiply(H1ar_y, temp_norm_h1_y[ndx])
    #  print(dotproduct_ar_normh1_ycomponent)
    H1dotproduct_ar_normh1_zcomponent = np.multiply(H1ar_z, temp_norm_h1_z[ndx])
    #  print(dotproduct_ar_normh1_zcomponent)
    H1dotproduct_ar_normh1_xycomponents = np.add(H1dotproduct_ar_normh1_xcomponent, H1dotproduct_ar_normh1_ycomponent)
    H1dotproduct_ar_normh1_allcomponents = np.add(H1dotproduct_ar_normh1_xycomponents, H1dotproduct_ar_normh1_zcomponent)
    #  print(dotproduct_ar_normh1_allcomponents)

    H1updated_theta_in_radians = np.zeros_like(H1theta_in_radians)

    if H1dotproduct_ar_normh1_allcomponents > 0:
        H1updated_theta_in_radians = H1theta_in_radians / 2.0
    else:
        H1updated_theta_in_radians = H1theta_in_radians / -2.0

    #  print("updated_theta_in_radians:", updated_theta_in_radians)

    H1updated_theta_in_degrees = np.degrees(H1updated_theta_in_radians)
    H1updated_theta_in_degrees_all.append(H1updated_theta_in_degrees.tolist())

    ###############################################################################

    #  Defining first quaternion based on the first axis of rotation (ar).
    H1q1_0 = np.cos(H1updated_theta_in_radians)
    #  print(q1_0)
    H1q1_1 = np.multiply(H1norm_ar_x, np.sin(H1updated_theta_in_radians))
    #  print(q1_1)
    H1q1_2 = np.multiply(H1norm_ar_y, np.sin(H1updated_theta_in_radians))
    #  print(q1_2)
    H1q1_3 = np.multiply(H1norm_ar_z, np.sin(H1updated_theta_in_radians))
    #  print(q1_3)

    e_0 = copy.deepcopy(H1q1_0)
    e_1 = copy.deepcopy(H1q1_1)
    e_2 = copy.deepcopy(H1q1_2)
    e_3 = copy.deepcopy(H1q1_3)

    ###############################################################################

    htemp_0 = np.multiply(np.subtract(np.add(np.square(e_0), np.square(e_1)), np.add(np.square(e_2), np.square(e_3))), temp_norm_h1_x[ndx])
    htemp_0 += np.multiply(2, np.multiply(np.add(np.multiply(e_1, e_2), np.multiply(e_0, e_3)), temp_norm_h1_y[ndx]))
    htemp_0 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_1, e_3), np.multiply(e_0, e_2)), temp_norm_h1_z[ndx]))
    #  print(htemp_0)

    htemp_1 = np.multiply(np.multiply(2, np.subtract(np.multiply(e_1, e_2), np.multiply(e_0, e_3))), temp_norm_h1_x[ndx])
    htemp_1 += np.multiply(np.add(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.subtract(np.multiply(e_2, e_2), np.multiply(e_3, e_3))), temp_norm_h1_y[ndx])
    htemp_1 += np.multiply(2, np.multiply(np.add(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), temp_norm_h1_z[ndx]))
    #  print(htemp_1)

    htemp_2 = np.multiply(np.multiply(2, np.add(np.multiply(e_1, e_3), np.multiply(e_0, e_2))), temp_norm_h1_x[ndx])
    htemp_2 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), temp_norm_h1_y[ndx]))
    htemp_2 += np.multiply(np.add(np.subtract(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.multiply(e_2, e_2)), np.multiply(e_3, e_3)), temp_norm_h1_z[ndx])
    #  print(htemp_2)

    h1_moved_x = copy.deepcopy(htemp_0)
    h1_moved_y = copy.deepcopy(htemp_1)
    h1_moved_z = copy.deepcopy(htemp_2)

    ###############################################################################

    htemp_0 = np.multiply(np.subtract(np.add(np.square(e_0), np.square(e_1)), np.add(np.square(e_2), np.square(e_3))), temp_norm_h2_x[ndx])
    htemp_0 += np.multiply(2, np.multiply(np.add(np.multiply(e_1, e_2), np.multiply(e_0, e_3)), temp_norm_h2_y[ndx]))
    htemp_0 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_1, e_3), np.multiply(e_0, e_2)), temp_norm_h2_z[ndx]))
    #  print(htemp_0)

    htemp_1 = np.multiply(np.multiply(2, np.subtract(np.multiply(e_1, e_2), np.multiply(e_0, e_3))), temp_norm_h2_x[ndx])
    htemp_1 += np.multiply(np.add(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.subtract(np.multiply(e_2, e_2), np.multiply(e_3, e_3))), temp_norm_h2_y[ndx])
    htemp_1 += np.multiply(2, np.multiply(np.add(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), temp_norm_h2_z[ndx]))
    #  print(htemp_1)

    htemp_2 = np.multiply(np.multiply(2, np.add(np.multiply(e_1, e_3), np.multiply(e_0, e_2))), temp_norm_h2_x[ndx])
    htemp_2 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), temp_norm_h2_y[ndx]))
    htemp_2 += np.multiply(np.add(np.subtract(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.multiply(e_2, e_2)), np.multiply(e_3, e_3)), temp_norm_h2_z[ndx])
    #  print(htemp_2)

    h2_after_h1_moved_x = copy.deepcopy(htemp_0)
    h2_after_h1_moved_y = copy.deepcopy(htemp_1)
    h2_after_h1_moved_z = copy.deepcopy(htemp_2)

    ###############################################################################

    H1z_mol_vect_x = np.subtract(np.multiply(h1_moved_y, h2_after_h1_moved_z), np.multiply(h1_moved_z, h2_after_h1_moved_y))
    H1z_mol_vect_y = np.subtract(np.multiply(h1_moved_z, h2_after_h1_moved_x), np.multiply(h1_moved_x, h2_after_h1_moved_z))
    H1z_mol_vect_z = np.subtract(np.multiply(h1_moved_x, h2_after_h1_moved_y), np.multiply(h1_moved_y, h2_after_h1_moved_x))

    H1z_mol_vect_x_squares = np.square(H1z_mol_vect_x)
    #  print(z_mol_vect_x_squares)
    H1z_mol_vect_y_squares = np.square(H1z_mol_vect_y)
    #  print(z_mol_vect_y_squares)
    H1z_mol_vect_z_squares = np.square(H1z_mol_vect_z)
    #  print(z_mol_vect_z_squares)
    H1z_mol_vect_xy_squares_added = np.add(H1z_mol_vect_x_squares, H1z_mol_vect_y_squares)
    H1z_mol_vect_all_squares_added = np.add(H1z_mol_vect_xy_squares_added, H1z_mol_vect_z_squares)
    #  print(z_mol_vect_all_squares_added)
    H1z_mol_vect_length = np.sqrt(H1z_mol_vect_all_squares_added)
    #  print(z_mol_vect_length)

    H1norm_z_mol_vect_x = np.divide(H1z_mol_vect_x, H1z_mol_vect_length)
    #  print(norm_z_mol_vect_x)
    H1norm_z_mol_vect_y = np.divide(H1z_mol_vect_y, H1z_mol_vect_length)
    #  print(norm_z_mol_vect_y)
    H1norm_z_mol_vect_z = np.divide(H1z_mol_vect_z, H1z_mol_vect_length)
    #  print(norm_z_mol_vect_z)

    ###############################################################################

    H1dotproduct_norm_z_mol_vect_zref_xcomponent = np.multiply(H1norm_z_mol_vect_x, z_ref[0])  # I know it's zero and that's because z_ref is [0,0,1]
    H1dotproduct_norm_z_mol_vect_zref_ycomponent = np.multiply(H1norm_z_mol_vect_y, z_ref[1])  # I know it's zero and that's because z_ref is [0,0,1]
    H1dotproduct_norm_z_mol_vect_zref_zcomponent = np.multiply(H1norm_z_mol_vect_z, z_ref[2])
    H1dotproduct_norm_z_mol_vect_zref_xycomponents = np.add(H1dotproduct_norm_z_mol_vect_zref_xcomponent, H1dotproduct_norm_z_mol_vect_zref_ycomponent)
    H1dotproduct_norm_z_mol_vect_zref_allcomponents = np.add(H1dotproduct_norm_z_mol_vect_zref_xycomponents, H1dotproduct_norm_z_mol_vect_zref_zcomponent)

    H1theta3p_in_radians = np.arccos(H1dotproduct_norm_z_mol_vect_zref_allcomponents)
    #  print(theta3p_in_radians)
    H1theta3p_in_degrees = np.degrees(H1theta3p_in_radians)
    #  print(theta3p_in_degrees)

    H1crossp_norm_z_mol_vect_z_ref_x = np.subtract(np.multiply(H1norm_z_mol_vect_y, z_ref[2]), np.multiply(H1norm_z_mol_vect_z, z_ref[1]))
    H1crossp_norm_z_mol_vect_z_ref_y = np.subtract(np.multiply(H1norm_z_mol_vect_z, z_ref[0]), np.multiply(H1norm_z_mol_vect_x, z_ref[2]))
    H1crossp_norm_z_mol_vect_z_ref_z = np.subtract(np.multiply(H1norm_z_mol_vect_x, z_ref[1]), np.multiply(H1norm_z_mol_vect_y, z_ref[0]))

    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent = np.multiply(H1crossp_norm_z_mol_vect_z_ref_x, h1_moved_x)
    #  print(dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent)
    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent = np.multiply(H1crossp_norm_z_mol_vect_z_ref_y, h1_moved_y)
    #  print(dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent)
    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent = np.multiply(H1crossp_norm_z_mol_vect_z_ref_z, h1_moved_z)
    #  print(dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent)
    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xycomponents = np.add(H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent, H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent)
    H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents = np.add(H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xycomponents, H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent)
    #  print(dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents)

    H1updated_theta3p_in_radians = np.zeros_like(H1theta3p_in_radians)

    if H1dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents < 0.0:
        H1updated_theta3p_in_radians = H1theta3p_in_radians / 2.0
    else:
        H1updated_theta3p_in_radians = H1theta3p_in_radians / -2.0

    #  print(updated_theta3p_in_radians)

    H1updated_theta3p_in_degrees = np.degrees(H1updated_theta3p_in_radians)
    H1updated_theta3p_in_degrees_all.append(H1updated_theta3p_in_degrees.tolist())
    ###############################################################################

    #  Defining second quaternion based on axis of rotation and angle.
    H1q2_0 = np.cos(H1updated_theta3p_in_radians)
    #  print(q2_0)
    H1q2_1 = np.multiply(x_ref[0], np.sin(H1updated_theta3p_in_radians))
    #  print(q2_1)
    H1q2_2 = np.multiply(x_ref[1], np.sin(H1updated_theta3p_in_radians))
    #  print(q2_2)
    H1q2_3 = np.multiply(x_ref[2], np.sin(H1updated_theta3p_in_radians))
    #  print(q2_3)

    ###############################################################################

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

    #  print(e_0)
    #  print(e_1)
    #  print(e_2)
    #  print(e_3)

    ##########################################################################

    H1final_q_0 = copy.deepcopy(e_0)
    H1final_q_1 = copy.deepcopy(e_1)
    H1final_q_2 = copy.deepcopy(e_2)
    H1final_q_3 = copy.deepcopy(e_3)

    H1final_q_0_all.append(H1final_q_0)
    H1final_q_1_all.append(H1final_q_1)
    H1final_q_2_all.append(H1final_q_2)
    H1final_q_3_all.append(H1final_q_3)

    #########################################################################################################################################################################################################################################################################################################################################################################################################
    # Transformation of H2
    H2ar_x = np.subtract(np.multiply(temp_norm_h2_y[ndx], x_ref[2]), np.multiply(temp_norm_h2_z[ndx], x_ref[1]))
    #  print("ar_x:", ar_x)
    H2ar_y = np.subtract(np.multiply(temp_norm_h2_z[ndx], x_ref[0]), np.multiply(temp_norm_h2_x[ndx], x_ref[2]))
    #  print("ar_y:", ar_y)
    H2ar_z = np.subtract(np.multiply(temp_norm_h2_x[ndx], x_ref[1]), np.multiply(temp_norm_h2_y[ndx], x_ref[0]))
    #  print("ar_z:", ar_z)

    #  Normalizing the first rotational axis (ar).

    H2ar_x_squares = np.square(H2ar_x)
    #  print(ar_x_squares)
    H2ar_y_squares = np.square(H2ar_y)
    #  print(ar_y_squares)
    H2ar_z_squares = np.square(H2ar_z)
    #  print(ar_z_squares)
    H2ar_xy_squares_added = np.add(H2ar_x_squares, H2ar_y_squares)
    H2ar_all_squares_added = np.add(H2ar_xy_squares_added, H2ar_z_squares)
    #  print(ar_all_squares_added)
    H2ar_length = np.sqrt(H2ar_all_squares_added)
    #  print("ar_length:", ar_length)

    H2norm_ar_x = np.nan_to_num(np.divide(H2ar_x, H2ar_length))
    #  print("norm_ar_x:", norm_ar_x)
    H2norm_ar_y = np.nan_to_num(np.divide(H2ar_y, H2ar_length))
    #  print("norm_ar_y:", norm_ar_y)
    H2norm_ar_z = np.nan_to_num(np.divide(H2ar_z, H2ar_length))
    #  print("norm_ar_z:", norm_ar_z)

    ###############################################################################

    #  Calculating theta from the dot product of x_ref vector and h1 vector
    #  https://math.stackexchange.com/questions/654315/how-to-convert-a-dot-product-of-two-vectors-to-the-angle-between-the-vectors

    H2dotproduct_h2_xref_xcomponent = np.multiply(temp_norm_h2_x[ndx], x_ref[0])
    H2dotproduct_h2_xref_ycomponent = np.multiply(temp_norm_h2_y[ndx], x_ref[1])  # I know it's zero and that's because x_ref is [1.0.0]
    H2dotproduct_h2_xref_zcomponent = np.multiply(temp_norm_h2_z[ndx], x_ref[2])  # I know it's zero and that's because x_ref is [1.0.0]
    H2dotproduct_h2_xref_xycomponents = np.add(H2dotproduct_h2_xref_xcomponent, H2dotproduct_h2_xref_ycomponent)
    H2dotproduct_h2_xref_allcomponents = np.add(H2dotproduct_h2_xref_xycomponents, H2dotproduct_h2_xref_zcomponent)
    #  print(dotproduct_h1_xref_allcomponents)
    H2theta_in_radians = np.arccos(H2dotproduct_h2_xref_allcomponents)
    #  print("theta_in_radians:", theta_in_radians)
    H2theta_in_degrees = np.degrees(H2theta_in_radians)
    #  print(theta_in_degrees)
    #  print(max(theta_in_degrees))
    #  print(min(theta_in_degrees))

    ###############################################################################

    #  Test to check for theta sign.

    H2dotproduct_ar_normh2_xcomponent = np.multiply(H2ar_x, temp_norm_h2_x[ndx])
    #  print(dotproduct_ar_normh1_xcomponent)
    H2dotproduct_ar_normh2_ycomponent = np.multiply(H2ar_y, temp_norm_h2_y[ndx])
    #  print(dotproduct_ar_normh1_ycomponent)
    H2dotproduct_ar_normh2_zcomponent = np.multiply(H2ar_z, temp_norm_h2_z[ndx])
    #  print(dotproduct_ar_normh1_zcomponent)
    H2dotproduct_ar_normh2_xycomponents = np.add(H2dotproduct_ar_normh2_xcomponent, H2dotproduct_ar_normh2_ycomponent)
    H2dotproduct_ar_normh2_allcomponents = np.add(H2dotproduct_ar_normh2_xycomponents, H2dotproduct_ar_normh2_zcomponent)
    #  print(dotproduct_ar_normh1_allcomponents)

    H2updated_theta_in_radians = np.zeros_like(H2theta_in_radians)

    if H2dotproduct_ar_normh2_allcomponents > 0:
        H2updated_theta_in_radians = H2theta_in_radians / 2.0
    else:
        H2updated_theta_in_radians = H2theta_in_radians / -2.0

    #  print("updated_theta_in_radians:", updated_theta_in_radians)

    H2updated_theta_in_degrees = np.degrees(H2updated_theta_in_radians)
    H2updated_theta_in_degrees_all.append(H2updated_theta_in_degrees.tolist())
    ###############################################################################

    #  Defining first quaternion based on the first axis of rotation (ar).
    H2q1_0 = np.cos(H2updated_theta_in_radians)
    #  print(q1_0)
    H2q1_1 = np.multiply(H2norm_ar_x, np.sin(H2updated_theta_in_radians))
    #  print(q1_1)
    H2q1_2 = np.multiply(H2norm_ar_y, np.sin(H2updated_theta_in_radians))
    #  print(q1_2)
    H2q1_3 = np.multiply(H2norm_ar_z, np.sin(H2updated_theta_in_radians))
    #  print(q1_3)

    e_0 = copy.deepcopy(H2q1_0)
    e_1 = copy.deepcopy(H2q1_1)
    e_2 = copy.deepcopy(H2q1_2)
    e_3 = copy.deepcopy(H2q1_3)

    ###############################################################################

    htemp_0 = np.multiply(np.subtract(np.add(np.square(e_0), np.square(e_1)), np.add(np.square(e_2), np.square(e_3))), temp_norm_h2_x[ndx])
    htemp_0 += np.multiply(2, np.multiply(np.add(np.multiply(e_1, e_2), np.multiply(e_0, e_3)), temp_norm_h2_y[ndx]))
    htemp_0 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_1, e_3), np.multiply(e_0, e_2)), temp_norm_h2_z[ndx]))
    #  print(htemp_0)

    htemp_1 = np.multiply(np.multiply(2, np.subtract(np.multiply(e_1, e_2), np.multiply(e_0, e_3))), temp_norm_h2_x[ndx])
    htemp_1 += np.multiply(np.add(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.subtract(np.multiply(e_2, e_2), np.multiply(e_3, e_3))), temp_norm_h2_y[ndx])
    htemp_1 += np.multiply(2, np.multiply(np.add(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), temp_norm_h2_z[ndx]))
    #  print(htemp_1)

    htemp_2 = np.multiply(np.multiply(2, np.add(np.multiply(e_1, e_3), np.multiply(e_0, e_2))), temp_norm_h2_x[ndx])
    htemp_2 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), temp_norm_h2_y[ndx]))
    htemp_2 += np.multiply(np.add(np.subtract(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.multiply(e_2, e_2)), np.multiply(e_3, e_3)), temp_norm_h2_z[ndx])
    #  print(htemp_2)

    h2_moved_x = copy.deepcopy(htemp_0)
    h2_moved_y = copy.deepcopy(htemp_1)
    h2_moved_z = copy.deepcopy(htemp_2)

    ###############################################################################

    htemp_0 = np.multiply(np.subtract(np.add(np.square(e_0), np.square(e_1)), np.add(np.square(e_2), np.square(e_3))), temp_norm_h1_x[ndx])
    htemp_0 += np.multiply(2, np.multiply(np.add(np.multiply(e_1, e_2), np.multiply(e_0, e_3)), temp_norm_h1_y[ndx]))
    htemp_0 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_1, e_3), np.multiply(e_0, e_2)), temp_norm_h1_z[ndx]))
    #  print(htemp_0)

    htemp_1 = np.multiply(np.multiply(2, np.subtract(np.multiply(e_1, e_2), np.multiply(e_0, e_3))), temp_norm_h1_x[ndx])
    htemp_1 += np.multiply(np.add(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.subtract(np.multiply(e_2, e_2), np.multiply(e_3, e_3))), temp_norm_h1_y[ndx])
    htemp_1 += np.multiply(2, np.multiply(np.add(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), temp_norm_h1_z[ndx]))
    #  print(htemp_1)

    htemp_2 = np.multiply(np.multiply(2, np.add(np.multiply(e_1, e_3), np.multiply(e_0, e_2))), temp_norm_h1_x[ndx])
    htemp_2 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), temp_norm_h1_y[ndx]))
    htemp_2 += np.multiply(np.add(np.subtract(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.multiply(e_2, e_2)), np.multiply(e_3, e_3)), temp_norm_h1_z[ndx])
    #  print(htemp_2)

    h1_after_h2_moved_x = copy.deepcopy(htemp_0)
    h1_after_h2_moved_y = copy.deepcopy(htemp_1)
    h1_after_h2_moved_z = copy.deepcopy(htemp_2)

    ###############################################################################

    H2z_mol_vect_x = np.subtract(np.multiply(h2_moved_y, h1_after_h2_moved_z), np.multiply(h2_moved_z, h1_after_h2_moved_y))
    H2z_mol_vect_y = np.subtract(np.multiply(h2_moved_z, h1_after_h2_moved_x), np.multiply(h2_moved_x, h1_after_h2_moved_z))
    H2z_mol_vect_z = np.subtract(np.multiply(h2_moved_x, h1_after_h2_moved_y), np.multiply(h2_moved_y, h1_after_h2_moved_x))

    H2z_mol_vect_x_squares = np.square(H2z_mol_vect_x)
    #  print(z_mol_vect_x_squares)
    H2z_mol_vect_y_squares = np.square(H2z_mol_vect_y)
    #  print(z_mol_vect_y_squares)
    H2z_mol_vect_z_squares = np.square(H2z_mol_vect_z)
    #  print(z_mol_vect_z_squares)
    H2z_mol_vect_xy_squares_added = np.add(H2z_mol_vect_x_squares, H2z_mol_vect_y_squares)
    H2z_mol_vect_all_squares_added = np.add(H2z_mol_vect_xy_squares_added, H2z_mol_vect_z_squares)
    #  print(z_mol_vect_all_squares_added)
    H2z_mol_vect_length = np.sqrt(H2z_mol_vect_all_squares_added)
    #  print(z_mol_vect_length)

    H2norm_z_mol_vect_x = np.divide(H2z_mol_vect_x, H2z_mol_vect_length)
    #  print(norm_z_mol_vect_x)
    H2norm_z_mol_vect_y = np.divide(H2z_mol_vect_y, H2z_mol_vect_length)
    #  print(norm_z_mol_vect_y)
    H2norm_z_mol_vect_z = np.divide(H2z_mol_vect_z, H2z_mol_vect_length)
    #  print(norm_z_mol_vect_z)

    ###############################################################################

    H2dotproduct_norm_z_mol_vect_zref_xcomponent = np.multiply(H2norm_z_mol_vect_x, z_ref[0])  # I know it's zero and that's because z_ref is [0,0,1]
    H2dotproduct_norm_z_mol_vect_zref_ycomponent = np.multiply(H2norm_z_mol_vect_y, z_ref[1])  # I know it's zero and that's because z_ref is [0,0,1]
    H2dotproduct_norm_z_mol_vect_zref_zcomponent = np.multiply(H2norm_z_mol_vect_z, z_ref[2])
    H2dotproduct_norm_z_mol_vect_zref_xycomponents = np.add(H2dotproduct_norm_z_mol_vect_zref_xcomponent, H2dotproduct_norm_z_mol_vect_zref_ycomponent)
    H2dotproduct_norm_z_mol_vect_zref_allcomponents = np.add(H2dotproduct_norm_z_mol_vect_zref_xycomponents, H2dotproduct_norm_z_mol_vect_zref_zcomponent)

    H2theta3p_in_radians = np.arccos(H2dotproduct_norm_z_mol_vect_zref_allcomponents)
    #  print(theta3p_in_radians)
    H2theta3p_in_degrees = np.degrees(H2theta3p_in_radians)
    #  print(theta3p_in_degrees)

    H2crossp_norm_z_mol_vect_z_ref_x = np.subtract(np.multiply(H2norm_z_mol_vect_y, z_ref[2]), np.multiply(H2norm_z_mol_vect_z, z_ref[1]))
    H2crossp_norm_z_mol_vect_z_ref_y = np.subtract(np.multiply(H2norm_z_mol_vect_z, z_ref[0]), np.multiply(H2norm_z_mol_vect_x, z_ref[2]))
    H2crossp_norm_z_mol_vect_z_ref_z = np.subtract(np.multiply(H2norm_z_mol_vect_x, z_ref[1]), np.multiply(H2norm_z_mol_vect_y, z_ref[0]))

    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent = np.multiply(H2crossp_norm_z_mol_vect_z_ref_x, h2_moved_x)
    #  print(dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent)
    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent = np.multiply(H2crossp_norm_z_mol_vect_z_ref_y, h2_moved_y)
    #  print(dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent)
    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent = np.multiply(H2crossp_norm_z_mol_vect_z_ref_z, h2_moved_z)
    #  print(dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent)
    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xycomponents = np.add(H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xcomponent, H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_ycomponent)
    H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents = np.add(H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_xycomponents, H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_zcomponent)
    #  print(dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents)

    H2updated_theta3p_in_radians = np.zeros_like(H2theta3p_in_radians)

    if H2dotproduct_crossp_norm_z_mol_vect_z_ref_h1_moved_allcomponents < 0.0:
        H2updated_theta3p_in_radians = H2theta3p_in_radians / 2.0
    else:
        H2updated_theta3p_in_radians = H2theta3p_in_radians / -2.0

    #  print(updated_theta3p_in_radians)

    H2updated_theta3p_in_degrees = np.degrees(H2updated_theta3p_in_radians)
    H2updated_theta3p_in_degrees_all.append(H2updated_theta3p_in_degrees.tolist())

    ###############################################################################

    #  Defining second quaternion based on axis of rotation and angle.
    H2q2_0 = np.cos(H2updated_theta3p_in_radians)
    #  print(q2_0)
    H2q2_1 = np.multiply(x_ref[0], np.sin(H2updated_theta3p_in_radians))
    #  print(q2_1)
    H2q2_2 = np.multiply(x_ref[1], np.sin(H2updated_theta3p_in_radians))
    #  print(q2_2)
    H2q2_3 = np.multiply(x_ref[2], np.sin(H2updated_theta3p_in_radians))
    #  print(q2_3)

    ###############################################################################

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

    #  print(e_0)
    #  print(e_1)
    #  print(e_2)
    #  print(e_3)

    ##########################################################################

    H2final_q_0 = copy.deepcopy(e_0)
    H2final_q_1 = copy.deepcopy(e_1)
    H2final_q_2 = copy.deepcopy(e_2)
    H2final_q_3 = copy.deepcopy(e_3)

    H2final_q_0_all.append(H2final_q_0)
    H2final_q_1_all.append(H2final_q_1)
    H2final_q_2_all.append(H2final_q_2)
    H2final_q_3_all.append(H2final_q_3)


#  TESTING MY ROTATIONS

e_0 = copy.deepcopy(H1final_q_0_all)
e_1 = copy.deepcopy(H1final_q_1_all)
e_2 = copy.deepcopy(H1final_q_2_all)
e_3 = copy.deepcopy(H1final_q_3_all)

#############################

htemp_0 = np.multiply(np.subtract(np.add(np.square(e_0), np.square(e_1)), np.add(np.square(e_2), np.square(e_3))), norm_h1_x)
htemp_0 += np.multiply(2, np.multiply(np.add(np.multiply(e_1, e_2), np.multiply(e_0, e_3)), norm_h1_y))
htemp_0 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_1, e_3), np.multiply(e_0, e_2)), norm_h1_z))
#  print(htemp_0)

htemp_1 = np.multiply(np.multiply(2, np.subtract(np.multiply(e_1, e_2), np.multiply(e_0, e_3))), norm_h1_x)
htemp_1 += np.multiply(np.add(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.subtract(np.multiply(e_2, e_2), np.multiply(e_3, e_3))), norm_h1_y)
htemp_1 += np.multiply(2, np.multiply(np.add(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), norm_h1_z))
#  print(htemp_1)

htemp_2 = np.multiply(np.multiply(2, np.add(np.multiply(e_1, e_3), np.multiply(e_0, e_2))), norm_h1_x)
htemp_2 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), norm_h1_y))
htemp_2 += np.multiply(np.add(np.subtract(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.multiply(e_2, e_2)), np.multiply(e_3, e_3)), norm_h1_z)
#  print(htemp_2)

FINAL_h1_moved_x = copy.deepcopy(htemp_0)
FINAL_h1_moved_y = copy.deepcopy(htemp_1)
FINAL_h1_moved_z = copy.deepcopy(htemp_2)

#############################

htemp_0 = np.multiply(np.subtract(np.add(np.square(e_0), np.square(e_1)), np.add(np.square(e_2), np.square(e_3))), norm_h2_x)
htemp_0 += np.multiply(2, np.multiply(np.add(np.multiply(e_1, e_2), np.multiply(e_0, e_3)), norm_h2_y))
htemp_0 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_1, e_3), np.multiply(e_0, e_2)), norm_h2_z))
#  print(htemp_0)

htemp_1 = np.multiply(np.multiply(2, np.subtract(np.multiply(e_1, e_2), np.multiply(e_0, e_3))), norm_h2_x)
htemp_1 += np.multiply(np.add(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.subtract(np.multiply(e_2, e_2), np.multiply(e_3, e_3))), norm_h2_y)
htemp_1 += np.multiply(2, np.multiply(np.add(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), norm_h2_z))
#  print(htemp_1)

htemp_2 = np.multiply(np.multiply(2, np.add(np.multiply(e_1, e_3), np.multiply(e_0, e_2))), norm_h2_x)
htemp_2 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), norm_h2_y))
htemp_2 += np.multiply(np.add(np.subtract(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.multiply(e_2, e_2)), np.multiply(e_3, e_3)), norm_h2_z)
#  print(htemp_2)

FINAL_h2_after_h1_moved_x = copy.deepcopy(htemp_0)
FINAL_h2_after_h1_moved_y = copy.deepcopy(htemp_1)
FINAL_h2_after_h1_moved_z = copy.deepcopy(htemp_2)

# print(max(FINAL_h1_moved_x))
# print(min(FINAL_h1_moved_x))
# print(max(FINAL_h1_moved_y))
# print(min(FINAL_h1_moved_y))
# print(max(FINAL_h1_moved_z))
# print(min(FINAL_h1_moved_z))
# print(max(FINAL_h2_after_h1_moved_x))
# print(min(FINAL_h2_after_h1_moved_x))
# print(max(FINAL_h2_after_h1_moved_y))
# print(min(FINAL_h2_after_h1_moved_y))
# print(max(FINAL_h2_after_h1_moved_z))
# print(min(FINAL_h2_after_h1_moved_z))

#############################################################################################################################################################################################################################################################################################
#  TESTING MY ROTATIONS

e_0 = copy.deepcopy(H2final_q_0_all)
e_1 = copy.deepcopy(H2final_q_1_all)
e_2 = copy.deepcopy(H2final_q_2_all)
e_3 = copy.deepcopy(H2final_q_3_all)

#############################

htemp_0 = np.multiply(np.subtract(np.add(np.square(e_0), np.square(e_1)), np.add(np.square(e_2), np.square(e_3))), norm_h2_x)
htemp_0 += np.multiply(2, np.multiply(np.add(np.multiply(e_1, e_2), np.multiply(e_0, e_3)), norm_h2_y))
htemp_0 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_1, e_3), np.multiply(e_0, e_2)), norm_h2_z))
#  print(htemp_0)

htemp_1 = np.multiply(np.multiply(2, np.subtract(np.multiply(e_1, e_2), np.multiply(e_0, e_3))), norm_h2_x)
htemp_1 += np.multiply(np.add(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.subtract(np.multiply(e_2, e_2), np.multiply(e_3, e_3))), norm_h2_y)
htemp_1 += np.multiply(2, np.multiply(np.add(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), norm_h2_z))
#  print(htemp_1)

htemp_2 = np.multiply(np.multiply(2, np.add(np.multiply(e_1, e_3), np.multiply(e_0, e_2))), norm_h2_x)
htemp_2 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), norm_h2_y))
htemp_2 += np.multiply(np.add(np.subtract(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.multiply(e_2, e_2)), np.multiply(e_3, e_3)), norm_h2_z)
#  print(htemp_2)

FINAL_h2_moved_x = copy.deepcopy(htemp_0)
FINAL_h2_moved_y = copy.deepcopy(htemp_1)
FINAL_h2_moved_z = copy.deepcopy(htemp_2)

#############################

htemp_0 = np.multiply(np.subtract(np.add(np.square(e_0), np.square(e_1)), np.add(np.square(e_2), np.square(e_3))), norm_h1_x)
htemp_0 += np.multiply(2, np.multiply(np.add(np.multiply(e_1, e_2), np.multiply(e_0, e_3)), norm_h1_y))
htemp_0 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_1, e_3), np.multiply(e_0, e_2)), norm_h1_z))
#  print(htemp_0)

htemp_1 = np.multiply(np.multiply(2, np.subtract(np.multiply(e_1, e_2), np.multiply(e_0, e_3))), norm_h1_x)
htemp_1 += np.multiply(np.add(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.subtract(np.multiply(e_2, e_2), np.multiply(e_3, e_3))), norm_h1_y)
htemp_1 += np.multiply(2, np.multiply(np.add(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), norm_h1_z))
#  print(htemp_1)

htemp_2 = np.multiply(np.multiply(2, np.add(np.multiply(e_1, e_3), np.multiply(e_0, e_2))), norm_h1_x)
htemp_2 += np.multiply(2, np.multiply(np.subtract(np.multiply(e_2, e_3), np.multiply(e_0, e_1)), norm_h1_y))
htemp_2 += np.multiply(np.add(np.subtract(np.subtract(np.multiply(e_0, e_0), np.multiply(e_1, e_1)), np.multiply(e_2, e_2)), np.multiply(e_3, e_3)), norm_h1_z)
#  print(htemp_2)

FINAL_h1_after_h2_moved_x = copy.deepcopy(htemp_0)
FINAL_h1_after_h2_moved_y = copy.deepcopy(htemp_1)
FINAL_h1_after_h2_moved_z = copy.deepcopy(htemp_2)

# print(max(FINAL_h2_moved_x))
# print(min(FINAL_h2_moved_x))
# print(max(FINAL_h2_moved_y))
# print(min(FINAL_h2_moved_y))
# print(max(FINAL_h2_moved_z))
# print(min(FINAL_h2_moved_z))
# print(max(FINAL_h1_after_h2_moved_x))
# print(min(FINAL_h1_after_h2_moved_x))
# print(max(FINAL_h1_after_h2_moved_y))
# print(min(FINAL_h1_after_h2_moved_y))
# print(max(FINAL_h1_after_h2_moved_z))
# print(min(FINAL_h1_after_h2_moved_z))

#############################################################################################################################################################################################################################################################################################

##########################################################################################################
# from IPython.display import clear_output  # to show in console the progress

# Initially the list to create a list of list with size N
N = len(H1final_q_0_all)
cutoff = 0.15000
Neighbors = [None] * N
for i in range(N):
    Neighbors[i] = []

for i in range(0, N - 1):
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

    # print("Current progress:", np.round(i / len(H1final_q_0_all) * 100, 2), "%")  # to show in console the progress


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



