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

import os,sys
import numpy as np
import random
sys.path.insert(0,"src/")
import reorder_jackknifes, filter_halos, sort_halos, jackknife_halos
import make_randoms, treecorr_interface, resum_cf
import build_delta_sigma

class catalog_analysis(object):
    """
    TODO: flow control, output paths for things,
    a set of mass bins, and a number of jackknife regions.
    """
    def __init__(self,dm_files,halo_file,cosmology,
                 treecorr_dict=None,flow_control=None,out_path="./",mass_bounds=None,
                 down_sampling=100):
        """
        dm_files: the base name for the LGADGET DM particle files.
        The files themselves have the extension .0 or .512 on them,
        corresponding to a jackknifed region (that has a random ordering).

        halo_file: the file that contains the halos identified by rockstar.

        treecorr_dict: a dictionary to be used to run treecorr in a customized way

        flow_control: used to leave out any parts of the process

        out_path: a path to the output directory

        mass_bounds: a custom split of halo masses
        """
        self.dm_files      = dm_files
        self.halo_file     = halo_file
        self.cosmology     = cosmology
        self.treecorr_dict = treecorr_dict
        self.flow_control  = flow_control
        self.out_path      = out_path
        self.mass_bounds   = mass_bounds
        self.down_sampling = down_sampling
        return

    def build_WL_signal(self):
        """
        This function calls all of the sub functions.
        """
        self.path_check()

        self.find_simulation_properties()

        self.mapping = reorder_jackknifes.reorder_jackknifes(self.ndm_jks,self.side_length,self.dm_files)

        #filter_halos.filter_halos(self.halo_file,self.filtered_halo_path)

        self.mass_bounds,self.log10_mass_bounds = sort_halos.sort_halos(self.filtered_halo_path,self.redshift,self.bounded_halo_path,self.mass_bounds)

        jackknife_halos.jackknife_halos(self.side_length,self.ndm_ndivs,self.redshift,self.log10_mass_bounds,self.bounded_halo_path,self.bounded_jk_output_directory_base,self.jk_halo_filename)

        #make_randoms.make_randoms(self.down_sampling,self.dm_count,self.ndm_jks,self.side_length,self.ndm_ndivs,self.rand_dm_path,self.rand_halo_path)

        #self.treecorr_dict = treecorr_interface.build_treecorr_dict(self.treecorr_dict,self.side_length,self.ndm_ndivs)

        os.system("mkdir -p %s"%self.cf_singles_jk_out_base%(self.redshift,self.log10_mass_bounds[0],self.log10_mass_bounds[1]))
        #treecorr_interface.run_treecorr(self.treecorr_dict,self.dm_files+".%d",self.bounded_jk_output_directory_base+self.jk_halo_filename,self.rand_dm_path,self.rand_halo_path,self.redshift,self.log10_mass_bounds,self.ndm_jks,self.ndm_ndivs,self.side_length,self.down_sampling,self.mapping,self.cf_singles_jk_out_base+self.cf_filename)

        os.system("mkdir -p %s"%self.cf_resum_out_base%(self.redshift,self.log10_mass_bounds[0],self.log10_mass_bounds[1]))
        os.system("mkdir -p %s"%self.cf_full_out_base%(self.redshift,self.log10_mass_bounds[0],self.log10_mass_bounds[1]))

        #resum_cf.resum_cf(self.cf_singles_jk_out_base+self.cf_filename,self.ndm_jks,self.ndm_ndivs,self.redshift,self.log10_mass_bounds,self.cf_resum_out_base+self.cf_resum_filename,self.cf_full_out_base+self.cf_full_filename)

        os.system("mkdir -p %s"%self.deltasigma_LOO_out_base%(self.redshift,self.log10_mass_bounds[0],self.log10_mass_bounds[1]))
        os.system("mkdir -p %s"%self.deltasigma_full_out_base%(self.redshift,self.log10_mass_bounds[0],self.log10_mass_bounds[1]))

        #GOT UP TO HERE
        build_delta_sigma.build_delta_sigma(self.cosmology,self.redshift,\
                                            self.log10_mass_bounds,\
                                            self.bounded_halo_path,\
                                            self.bounded_jk_output_directory_base+self.jk_halo_filename,\
                                            self.cf_full_out_base+self.cf_full_filename,\
                                            self.cf_resum_out_base+self.cf_resum_filename,\
                                            self.deltasigma_LOO_out_base+self.deltasigma_LOO_filename,\
                                            self.deltasigma_full_out_base+self.deltasigma_full_filename)

        #self.build_delta_sigmas()
        #self.resum_delta_sigmas()

    def path_check(self):
        """
        This checks to see if the paths to the dark matter
        and halo data exist as well as creates
        the output path.
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

        print "Creating output directory."
        self.bounded_halo_path = self.out_path+"halos_z%.2f_%.1f_%.1f.txt"
        self.jk_out_path = self.out_path+"/jackknifed_halos/"
        self.bounded_jk_output_directory_base = self.jk_out_path+"jackknifed_halos_z%.2f_%.1f_%.1f/"
        self.jk_halo_filename = "halos_z%.2f_%.1f_%.1f_jk%d.txt"
        self.out_rand_path = self.out_path+"/randoms/"
        os.system("mkdir -p %s"%self.out_path)
        os.system("mkdir -p %s"%self.jk_out_path)
        os.system("mkdir -p %s"%self.out_rand_path)

        halo_filename = self.halo_file.split("/")[-1]
        self.filtered_halo_path = self.out_path+"filtered_%s.txt"%halo_filename
        DS = self.down_sampling
        out_dm_name,out_halo_name = "/rand_dm_DS%d.txt"%DS,"/rand_halo_DS%d.txt"%DS
        self.rand_dm_path, self.rand_halo_path = self.out_rand_path+out_dm_name,self.out_rand_path+out_halo_name
        
        self.cf_jk_out_base = self.out_path+"/jackknifed_CF/"
        os.system("mkdir -p %s"%self.cf_jk_out_base)
        self.cf_singles_jk_out_base = self.cf_jk_out_base+"/cf_z%.2f_%.1f_%.1f_singles/"
        self.cf_filename = "/cf_z%.2f_%.1f_%.1f_dmjk%d_halojk%d.txt"
        self.cf_resum_out_base = self.cf_jk_out_base+"/cf_z%.2f_%.1f_%.1f_resum/"
        self.cf_resum_filename = "/cf_resum_z%.2f_%.1f_%.1f_jk%d.txt"
        self.cf_full_out_base  = self.cf_jk_out_base+"/cf_z%.2f_%.1f_%.1f_full/"
        self.cf_full_filename  = "/cf_full_z%.2f_%.1f_%.1f.txt"

        self.deltasigma_jk_out_base = self.out_path+"/jackknifed_DS/"
        os.system("mkdir -p %s"%(self.deltasigma_jk_out_base))
        self.deltasigma_LOO_out_base = self.deltasigma_jk_out_base+"/DS_z%.2f_%.1f_%.1f_LOO/"
        self.deltasigma_LOO_filename = self.deltasigma_jk_out_base+"/DS_LOO_z%.2f_%.1f_%.1f.txt"
        self.deltasigma_full_out_base = self.deltasigma_jk_out_base+"/DS_z%.2f_%.1f_%.1f_full/"
        self.deltasigma_full_filename = self.deltasigma_jk_out_base+"/DS_full_z%.2f_%.1f_%.1f.txt"

        print "\tOutput directories created"
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
        self.ndm_ndivs = int(round(self.ndm_jks**(1./3.)))
        self.redshift = header['redshift']
        self.hubble_const = header['h']
        print "\tSnapshot header read. Simulation properties saved."
        return

if __name__ == '__main__':
    cosmo = {"h":0.7,"om":0.3,"ok":0.0,"ode":0.7,"sigma8":0.8,"ns":0.96}

    test = catalog_analysis("test_data/dm_files/snapshot_000",\
                            "test_data/halo_files/outbgc2_0.list",\
                            cosmology=cosmo,out_path="./output/",down_sampling=10)
    test.build_WL_signal()

