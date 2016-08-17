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
        self.dm_files       = dm_files
        self.halo_file     = halo_file
        self.treecorr_dict = treecorr_dict
        self.flow_control  = flow_control
        self.out_path     = out_path
        return

    def build_WL_signal(self):
        """
        This function calls all of the sub functions.
        """
        self.path_check()
        self.find_simulation_properties()
        #self.reorder_jackknifes()
        #self.sort_halos()
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
        #Check that the files exist
        if not os.path.exists(self.dm_files):
            raise Exception("Error: dark matter files don't exst.")
        if not os.path.exists(self.halo_file):
            raise Exception("Error: halo files don't exst.")
        print "\tDM and halo paths exist."
        return

    def find_simulation_properties(self):
        import pygadgetreader as pgr
        header = pgr.readheader(self.dm_files,"header")
        self.pgrheader = header
        self.particle_mass = max(header['massTable'])
        self.side_length = header['boxsize']
        self.dm_count = header['ndm']
        self.ndm_jks = header['nfiles']
        self.redshift = header['redshift']
        self.hubble_const = header['h']
        print header
        return

if __name__ == '__main__':
    test = catalog_analysis("test_data/","test_data/")
    test.build_WL_signal()
