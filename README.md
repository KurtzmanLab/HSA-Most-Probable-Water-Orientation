# HSA-Most-Probable-Water-Orientation
Code to find the most probable orientations in each HSA cluster

The following two python scripts need to be run while inside the produced SSTMap_HSA directory


water_orientation_clustering.py is a python script that takes each HSA cluster file such as 'cluster.000001.pdb' and find the most probable water orientations in that cluster.
The script outputs all most probable water orientations in a pdb file named as 'Water_orientational_clusters_of_cluster.000001.pdb', for example.
The script adds a popoulation percentage of each water orientation next to the coordinates of the oxygen atom.

Example for usage is 
$water_orientation_clustering.py cluster.000001.pdb

orientational_cluster_centering.py is a python script that takes the pdb file that contains all most probable waters and centers these waters to the position of the cluster oxygen found in 'clustercenterfile.pdb' and outputs them in a file named 'Centered_Water_orientational_clusters_of_cluster.000001.pdb'

Example for usage is
$orientational_cluster_centering.py Water_orientational_clusters_of_cluster.000001.pdb


