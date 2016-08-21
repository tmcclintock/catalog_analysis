"""
This file contains stand-alone versions of the functions
needed to interface properly with TreeCorr,
including the build_treecorr_dict()
function and the run_treecorr() function.
"""
import treecorr
import numpy as np

def build_treecorr_dict(treecorr_dict,side_length,ndivs):
    if treecorr_dict is not None:
        print "Using premade TreeCorr dictionary."
        return treecorr_dict
    print "Building TreeCorr dictionary."
    config = {}
    step = side_length/ndivs
    config['nbins'] = 50 #arbitrary
    config['min_sep'] = 0.1 # Mpc/h presumably; arbitrary
    config['max_sep'] = 50.0 # Mpc/h presumably; arbitrary
    print "\tTreeCorr dictionary built with nbins = %d, min_sep = %.1e, max_sep = %.1e"%(config['nbins'],config['min_sep'],config['max_sep'])
    return config

def run_treecorr(config,dm_base_file_path,halo_base_file_path,dm_random_base,halo_random_base,redshift,mass_bounds,N_jks,ndivs):
    return
