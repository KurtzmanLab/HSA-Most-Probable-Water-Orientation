#!/usr/bin/python

# ===============================================================================
#
#          FILE: orientational_cluster_centering.py
#
#         USAGE: ./orientational_cluster_centering.py 'name_of_cluster.pdb'
#
#   DESCRIPTION:
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Mossa Ghattas, mossa.ghattas@gmail.com
#  ORGANIZATION: Lehman College, Chemistry Department
#       CREATED: 04/03/2020 18:10:42
#      REVISION:  ---
# ====================================================================


import math
import numpy as np
import sys
import copy
import datetime
import os

centerfile = "clustercenterfile.pdb"
inputfile = sys.argv[1]
hsa_num = int(inputfile[40:46]) + 1

with open(centerfile, "r") as file:  # opens the file and reads each line as an item in a list and closes it.

    np.seterr(divide='ignore', invalid='ignore')  # https://stackoverflow.com/questions/14861891/runtimewarning-invalid-value-encountered-in-divide

    in_file = file.readlines()  # converts what is read to a list

    oxygen = in_file[hsa_num]  # extracting oxygen from in_file

center_np_o_x = float(oxygen[32:38])
center_np_o_y = float(oxygen[40:46])
center_np_o_z = float(oxygen[48:54])


with open(inputfile, "r") as file:  # opens the file and reads each line as an item in a list and closes it.

    np.seterr(divide='ignore', invalid='ignore')  # https://stackoverflow.com/questions/14861891/runtimewarning-invalid-value-encountered-in-divide

    in_file = file.readlines()  # converts what is read to a list

    end = len(in_file)

    oxygen = in_file[2:end:3]  # extracting oxygen from in_file
    h1 = in_file[3:end:3]  # extracting h1 from in_file
    h2 = in_file[4:end:3]  # extracting h2 from in_file

oxygen_x = []
for i in oxygen:
    oxygen_x.append(float(i[32:38]))
np_o_x = np.array(oxygen_x)


oxygen_y = []
for j in oxygen:
    oxygen_y.append(float(j[40:46]))
np_o_y = np.array(oxygen_y)


oxygen_z = []
for k in oxygen:
    oxygen_z.append(float(k[48:54]))
np_o_z = np.array(oxygen_z)


percentage_of_clusters = []
for z in oxygen:
    percentage_of_clusters.append(str(round(float(z[56:59]), 1)))


h1_x = []
for i in h1:
    h1_x.append(float(i[32:38]))
np_h1_x = np.array(h1_x)


h1_y = []
for j in h1:
    h1_y.append(float(j[40:46]))
np_h1_y = np.array(h1_y)


h1_z = []
for k in h1:
    h1_z.append(float(k[48:54]))
np_h1_z = np.array(h1_z)


h2_x = []
for i in h2:
    h2_x.append(float(i[32:38]))
np_h2_x = np.array(h2_x)


h2_y = []
for j in h2:
    h2_y.append(float(j[40:46]))
np_h2_y = np.array(h2_y)


h2_z = []
for k in h2:
    h2_z.append(float(k[48:54]))
np_h2_z = np.array(h2_z)

###############################################################################

diff_cluster_individual_x = np.subtract(np_o_x, center_np_o_x)
diff_cluster_individual_y = np.subtract(np_o_y, center_np_o_y)
diff_cluster_individual_z = np.subtract(np_o_z, center_np_o_z)

individual_fixed_to_center_o_x = np.subtract(np_o_x, diff_cluster_individual_x)
individual_fixed_to_center_o_y = np.subtract(np_o_y, diff_cluster_individual_y)
individual_fixed_to_center_o_z = np.subtract(np_o_z, diff_cluster_individual_z)

individual_fixed_to_center_h1_x = np.subtract(np_h1_x, diff_cluster_individual_x)
individual_fixed_to_center_h1_y = np.subtract(np_h1_y, diff_cluster_individual_y)
individual_fixed_to_center_h1_z = np.subtract(np_h1_z, diff_cluster_individual_z)

individual_fixed_to_center_h2_x = np.subtract(np_h2_x, diff_cluster_individual_x)
individual_fixed_to_center_h2_y = np.subtract(np_h2_y, diff_cluster_individual_y)
individual_fixed_to_center_h2_z = np.subtract(np_h2_z, diff_cluster_individual_z)
###############################################################################

###############################################################################

outfile = open("Centered_{}".format(inputfile), "w")
outfile.write("REMARKS\n")
# outfile.write("This pdb file has the most probable orientations in a water cluster.\n")
outfile.write("The centered orientational cluster.\n")
atmndx = 1
for i in np.arange(len(individual_fixed_to_center_o_x)):
    outfile.write("{:6s}{:5d} {:^4s} WT{:1s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}  {}%\n".format("ATOM", atmndx, "OW", str(i+1), "O", hsa_num-1, individual_fixed_to_center_o_x[i], individual_fixed_to_center_o_y[i], individual_fixed_to_center_o_z[i], percentage_of_clusters[i]))
    atmndx += 1
    outfile.write("{:6s}{:5d} {:^4s} WT{:1s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format("ATOM", atmndx, "H1", str(i+1), "H", hsa_num-1, individual_fixed_to_center_h1_x[i], individual_fixed_to_center_h1_y[i], individual_fixed_to_center_h1_z[i]))
    atmndx += 1
    outfile.write("{:6s}{:5d} {:^4s} WT{:1s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format("ATOM", atmndx, "H2", str(i+1), "H", hsa_num-1, individual_fixed_to_center_h2_x[i], individual_fixed_to_center_h2_y[i], individual_fixed_to_center_h2_z[i]))
    atmndx += 1
outfile.close()

###############################################################################

