"""
This is the top level object that contains the interface to the 
catalog analysis code. This code recieves a rockstar simulation
and performs certain calculations on it.

For more information visit:
https://github.com/tmcclintock/catalog_analysis

Required dependencies:
pygadgetreader
treecorr
numpy
scipy
"""

import os
import numpy as np

class catalog_analysis(object):
    """
    TODO: flow control, output paths for things,
    a set of mass bins, and a number of jackknife regions.
    """
    def __init__(self,dm_files,halo_file,
                 treecorr_dict=None,flow_control=None,out_path="./"):
        """
        dm_files: the base name for the LGADGET DM particle files.
        The files themselves have the extension .0 or .512 on them,
        corresponding to a jackknifed region (that has a random ordering).

        halo_file: the file that contains the halos identified by rockstar.
        """
        self.dm_files      = dm_files
        self.halo_file     = halo_file
        self.treecorr_dict = treecorr_dict
        self.flow_control  = flow_control
        self.out_path      = out_path
        return

    def build_WL_signal(self):
        """
        This function calls all of the sub functions.
        """
        self.path_check()
        self.find_simulation_properties()
        self.reorder_jackknifes()
        #self.sort_halos()
        #self.make_randoms()
        #self.run_treecorr()
        #self.resum_correlation_functions()
        #self.build_delta_sigmas()
        #self.resum_delta_sigmas()

    def path_check(self):
        """
        This checks to see if the paths to the dark matter
        and halo data exist.
        """
        print "Checking DM and halo paths."
        parts = self.dm_files.split("/")
        self.dm_filename = parts[-1]
        if len(parts)==1:
            self.dm_path = "./"
        else:
            self.dm_path = parts[0]
            for p in parts[1:-1]:
                self.dm_path = self.dm_path+"/"+p
        #Check that the files exist
        if not os.path.exists(self.dm_path):
            raise Exception("dark matter path don't exst.")
        if not os.path.exists(self.halo_file):
            raise Exception("halo files don't exst.")
        print "\tDM and halo paths exist."
        return

    def find_simulation_properties(self):
        """
        This finds the attributes of the simulation such as
        side length, number of DM paticles, number of DM JK files,
        the redshift, the hubble constant.
        TODO: figure out the cosmology.
        """
        print "Reading the snapshot header."
        import pygadgetreader as pgr
        header = pgr.readheader(self.dm_files,"header")
        self.pgrheader = header
        self.particle_mass = max(header['massTable'])
        self.side_length = header['boxsize']
        self.dm_count = header['ndm']
        self.ndm_jks = header['nfiles']
        self.redshift = header['redshift']
        self.hubble_const = header['h']
        print "\tSnapshot header read. Simulation properties saved."
        return

    def reorder_jackknifes(self):
        """
        This reads in the snapshot DM JK files and
        creates a dictionary that pairs the LGADGET JK number
        with our ordering.
        """
        print "Finding the mapping of LGADGET jackknifes to our jackknifes."
        import pygadgetreader as pgr
        N = self.ndm_jks #Number of LGADGET JK files
        ndivs = int(round(N**(1./3.)))
        self.ndm_ndivs = ndivs
        side = self.side_length
        step = side/ndivs
        mapping = {}
        for i in xrange(0,2): #N):
            pos = pgr.readsnap(self.dm_files+".%d"%i,"pos","dm",single=True,suppress=True)
            xi,yi,zi = [int(q) for q in np.mean(pos,0)/step]
            jkindex = zi + ndivs*yi + ndivs*ndivs*xi
            mapping[jkindex] = i
        self.mapping = mapping
        print "\tMapping complete on %d files."%N
        return

    def sort_halos(self):
        out_path = self.out_path
        

if __name__ == '__main__':
    test = catalog_analysis("test_data/dm_files/snapshot_000","test_data/halo_files/outbgc2_0.list")
    test.build_WL_signal()
