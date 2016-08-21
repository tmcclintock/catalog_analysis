"""
This is a stand-alone version of the make_randoms function that
appeared originally in the catalog_analysis.py file.
"""
import numpy as np

def make_randoms(DS,N_parts,N_jks,side_length,ndivs,rand_dm_path,rand_halo_path):
    print "Creating random points."
    step = side_length/ndivs
    N_dm_jkparts = N_parts/N_jks/DS
    N_halo_jkparts = N_dm_jkparts/2
    rand_dm = np.random.random((N_dm_jkparts,3))*step
    rand_halos = np.random.random((N_halo_jkparts,3))*step
    np.savetxt(rand_dm_path,rand_dm)
    np.savetxt(rand_halo_path,rand_halos)
    print "\tRandom points created."
