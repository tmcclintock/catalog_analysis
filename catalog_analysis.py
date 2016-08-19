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
                 treecorr_dict=None,flow_control=None,out_path="./",mass_bounds=None):
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
        self.treecorr_dict = treecorr_dict
        self.flow_control  = flow_control
        self.out_path      = out_path
        self.mass_bounds   = mass_bounds
        return

    def build_WL_signal(self):
        """
        This function calls all of the sub functions.
        """
        self.path_check()
        self.find_simulation_properties()
        self.reorder_jackknifes()
        #self.filter_halos()
        #self.sort_halos()
        #self.jackknife_halos()
        self.make_randoms()
        #self.run_treecorr()
        #self.resum_correlation_functions()
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
        out_path = self.out_path
        self.jk_out_path = self.out_path+"/jackknifed_halos/"
        os.system("mkdir -p %s"%out_path)
        os.system("mkdir -p %s"%self.jk_out_path)
        halo_filename = self.halo_file.split("/")[-1]
        self.filtered_halo_path = self.out_path+"filtered_%s.txt"%halo_filename
        print "\tOutput directory created"
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

    def filter_halos(self):
        """
        This function takes in the halo catalog and filters out the halos
        that are either subhalos (pid > -1.0) or aren't massive
        enouch (Np >= 200).
        """
        print "Filtering halo list."
        data = open(self.halo_file,"r")
        header = data.readline()
        outdata = open(self.filtered_halo_path,"w")
        outdata.write("#M X Y Z\n")
        for line in data:
            ID,DID,M,Vmax,Vrms,R200,Rs,Np,x,y,z,vx,vy,vz,pid = [float(item) for item in line.split()]
            if pid < 0.0 and Np >= 200:
                outdata.write("%e %e %e %e\n"%(M,x,y,z))
        outdata.close()
        print "\tHalos filtered."
        return

    def sort_halos(self):
        """
        This function takes the filtered halo file and
        sorts the halos into mass bins.
        """
        print "Sorting filtered halos."
        halos = np.genfromtxt(self.filtered_halo_path)
        Mass,x,y,z = halos.T
        lM = np.log10(Mass)
        if self.mass_bounds is None:
            Mmin,Mmax,Nbins = min(Mass),max(Mass),3
            edges = np.linspace(np.log10(Mmin),np.log10(Mmax),Nbins+1)
            self.mass_bounds = np.array([edges[:-1],edges[1:]]).T
            np.savetxt(self.out_path+"mass_bound_list.txt",self.mass_bounds)
            print "\tCreating %d mass bins between %.2e and %.2e"%(Nbins,Mmin,Mmax)
        outlist = []
        for i in range(len(self.mass_bounds)):
            mbs = self.mass_bounds[i]
            mbs_path = self.out_path+"halos_z%.2f_MB%d.txt"%(self.redshift,i)
            outlist.append(open(mbs_path,"w"))

        for i in range(len(halos)):
            for j in range(len(self.mass_bounds)):
                mbs = self.mass_bounds[j]
                if lM[i] >= mbs[0] and lM[i]< mbs[1]:
                    outlist[j].write("%e %e %e %e\n"%(Mass[i],x[i],y[i],z[i]))
                    break
                continue
            if lM[i] == self.mass_bounds[-1,1]: #The edge case
                outlist[-1].write("%e %e %e %e\n"%(Mass[i],x[i],y[i],z[i]))
        for i in range(len(self.mass_bounds)):
            outlist[i].close()
        print "\tHalos are now sorted into mass bins."
        return

    def jackknife_halos(self):
        """
        This function takes the sorted halos and
        jackknifes them based on how many divisions there are
        in the DM files.
        """
        print "Jackknifing halo files."
        side = self.side_length
        ndivs = self.ndm_ndivs
        step = side/ndivs
        outpath = self.jk_out_path
        for i in range(len(self.mass_bounds)):
            mbs_dirname = "jackknifed_halos_z%.2f_MB%d/"%(self.redshift,i)
            os.system("mkdir -p %s"%(outpath+mbs_dirname))
            
            #Read in the halos
            halos = np.genfromtxt(self.out_path+"halos_z%.2f_MB%d.txt"%(self.redshift,i))
            M,x,y,z = halos.T
            xi = np.floor(x/step)
            yi = np.floor(y/step)
            zi = np.floor(z/step)
            jkindices = zi + ndivs*yi + ndivs*ndivs*xi

            #Open all of the jackknife files
            outlist = []
            for jk in range(self.ndm_jks):
                outlist.append(open(outpath+mbs_dirname+"halos_z%.2f_MB%d_jk%d.txt"%(self.redshift,i,jk),"w"))
            
            
            for jki in range(len(jkindices)):
                outlist[int(jkindices[jki])].write("%e %e %e %e\n"%(M[jki],x[jki],y[jki],z[jki]))            
            for jk in range(self.ndm_jks):
                outlist[jk].close()
            continue #end i in mass_bounds
        print "\tHalos are now jackknifed."
        return

    def make_randoms(self):
        """
        This function creates the random points used in the
        correlation function pair counting calculation.
        """
        print "Creating random points."
        
        print "\tRandom points created."
        return
if __name__ == '__main__':
    test = catalog_analysis("test_data/dm_files/snapshot_000","test_data/halo_files/outbgc2_0.list",out_path="./output/")
    test.build_WL_signal()
