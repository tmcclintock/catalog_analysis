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
import reorder_jackknifes, filter_halos

class catalog_analysis(object):
    """
    TODO: flow control, output paths for things,
    a set of mass bins, and a number of jackknife regions.
    """
    def __init__(self,dm_files,halo_file,
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
        filter_halos.filter_halos(self.halo_file,self.filtered_halo_path)
        self.sort_halos()
        #self.jackknife_halos()
        #self.make_randoms()
        #if self.treecorr_dict is None:
        #    self.build_treecorr_dict()
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
        mbs_base = "jackknifed_halos_z%.2f_MB%d/"
        self.mbs_base = mbs_base
        jk_halo_fname = self.jk_out_path+mbs_base+"halos_z%.2f_MB%d_jk%d.txt"
        self.jk_halo_fname = jk_halo_fname

        out_rand = self.out_path+"/randoms/"
        self.out_rand = out_rand
        os.system("mkdir -p %s"%out_path)
        os.system("mkdir -p %s"%self.jk_out_path)
        os.system("mkdir -p %s"%out_rand)
        halo_filename = self.halo_file.split("/")[-1]
        self.filtered_halo_path = self.out_path+"filtered_%s.txt"%halo_filename
        DS = self.down_sampling
        out_dm_name,out_halo_name = "/rand_dm_DS%d.txt"%DS,"/rand_halo_DS%d.txt"%DS
        self.rand_dm_path, self.rand_halo_path = out_rand+out_dm_name,out_rand+out_halo_name
        
        cf_jk_out_base = out_path+"/jackknifed_CF/"
        self.cf_jk_out_base = cf_jk_out_base
        os.system("mkdir -p %s"%cf_jk_out_base)
        cf_MB_singles_jk_out_base = cf_jk_out_base+"/cf_z%.2f_MB%d_singles/"
        self.cf_MB_singles_jk_out_base = cf_MB_singles_jk_out_base
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
        self.ndm_ndivs = int(round(self.ndm_jks**(1./3.)))
        self.redshift = header['redshift']
        self.hubble_const = header['h']
        print "\tSnapshot header read. Simulation properties saved."
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
        mbs_base = self.mbs_base
        halo_jk_fname = self.jk_halo_fname
        redshift = self.redshift
        for i in range(len(self.mass_bounds)):
            mbs_dirname = mbs_base%(self.redshift,i)
            os.system("mkdir -p %s"%(outpath+mbs_dirname))
            
            #Read in the halos
            halos = np.genfromtxt(self.out_path+"halos_z%.2f_MB%d.txt"%(redshift,i))
            M,x,y,z = halos.T
            xi = np.floor(x/step)
            yi = np.floor(y/step)
            zi = np.floor(z/step)
            jkindices = zi + ndivs*yi + ndivs*ndivs*xi

            #Open all of the jackknife files
            outlist = []
            for jk in range(self.ndm_jks):
                outlist.append(open(halo_jk_fname%(redshift,i,redshift,i,self.jk),"w"))
            
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
        It is during the treecorr call stage that the individual
        random segments are translated to the correct jackknife
        regions
        """
        print "Creating random points."
        DS = self.down_sampling
        out_dm_name, out_halo_name = self.rand_dm_path, self.rand_halo_path
        Nparts = self.dm_count
        N_jks = self.ndm_jks
        side = self.side_length
        ndivs = self.ndm_ndivs
        step = side/ndivs
        N_dm_jkparts = Nparts/N_jks/DS
        N_halo_jkparts = N_dm_jkparts/2
        x,y,z = np.random.random((3,N_dm_jkparts))*step
        rand_dm = np.array([x,y,z]).T
        xh,yh,zh = np.random.random((3,N_halo_jkparts))*step
        rand_halos = np.array([xh,yh,zh]).T
        np.savetxt(out_dm_name,rand_dm)
        np.savetxt(out_halo_name,rand_halos)
        print "\tRandom points created."
        return

    def build_treecorr_dict():
        """
        This builds our treecorr dictionary in preparation
        for running treecorr.
        """
        print "Creating treecorr dictionary."
        N_jks, side, ndivs = self.ndm_jks, self.side_length, self.ndm_ndivs
        step = side/ndivs
        config = {}
        config['nbins'] = 50 # arbitrary
        config['min_sep'] = 0.1 # Mpc/h
        config['max_sep'] = 50.0# Mpc/h
        

        print "\tTreecorr dictionary built."
        return 


    def run_treecorr(self):
        """
        This function takes the DM snapshots and random files and loops 
        over halo subsets to calculate the halo-matter correlation
        function, xi(r), as well as the rest of the quantities that come with it.

        Note: this creates the single correlation/cross-correlation files.
        It does not perform the re-summation to get the leave-one-out (LOO)
        CFs.
        """
        print "Reading in all randoms, dark matter and halo files."
        import pygadgetreader as pgr
        import treecorr
        redshift = self.redshift
        rand_dm_path, rand_halo_path = self.rand_dm_path, self.rand_halo_path
        N_jks, side, ndivs = self.ndm_jks, self.side_length, self.ndm_ndivs
        step = side/ndivs
        mapping = self.mapping
        redshift = self.redshift
        jk_halo_fname = self.jk_halo_fname
        dm_path = self.dm_files+".%d"
        cf_jk_out_base = self.cf_jk_out_base
        cf_MB_singles_jk_out_base = self.cf_MB_singles_jk_out_base
        mass_bounds = self.mass_bounds

        #Read in the randoms
        rand_dm = np.loadtxt(rand_dm_path)
        rand_halo = np.loadtxt(rand_halo_path)
        print "\tRandom DM shape:",rand_dm.shape
        print "\tRandom halo shape:",rand_halo.shape

        #Loop over all JK pairs, first loop over the DM indices
        for i in xrange(0,1):#ndivs):
            for j in xrange(0,1):#ndivs):
                for k in xrange(0,1):#ndivs):
                    unmapped_dm_jkindex = k + ndivs*j + ndivs*ndivs*i
                    dm_jkindex = mapping[unmapped_dm_jkindex]
                    dm_parts_all = pgr.readsnap(dm_path%dm_jkindex,"pos","dm",single=True,suppress=True)
                    DS, Ndm = self.down_sampling, len(dm_parts_all)
                    fkeep = 1./DS
                    dm_parts = []
                    for dmp in range(len(dm_parts_all)):
                        r = np.random.random()
                    if r <= fkeep:
                        dm_parts.append(dm_parts_all[dmp])
                    dm_parts = np.array(dm_parts)
                    print "\tDM particles jk%d all/DS:"%dm_jkindex, dm_parts_all.shape, dm_parts.shape

                    #Move the DM randoms
                    rand_dm[:,0] += i*step #x
                    rand_dm[:,1] += j*step #y
                    rand_dm[:,2] += k*step #z

                    #Now loop over halo indices
                    for ii in xrange(0,1):#ndivs):
                        for jj in xrange(0,1):#ndivs):
                            for kk in xrange(0,1):#ndivs):
                                halo_jkindex = kk + ndivs*jj + ndivs*ndivs*ii
                                
                                #Translate the halo randoms over
                                rand_halo[:,0] += ii*step #x
                                rand_halo[:,1] += jj*step #y
                                rand_halo[:,2] += kk*step #z

                                for mbs_index in range(0,1):#len(mass_bounds)):
                                    halos = np.atleast_2d(np.loadtxt(jk_halo_fname%(redshift,mbs_index,redshift,mbs_index,halo_jkindex)))
                                    halo_pos = halos[:,1:]
                                    print "\tHalos MB%d jk%d shape:"%(mbs_index,halo_jkindex),halo_pos.shape

                                    #Set up treecorr catalog
                                    print "\tInterfacing with TreeCorr to calculate xi(r)."
                                        #dm_cat = treecorr.Catalog(x=dm_parts[:,0],y=dm_parts[:,1],z=dm_parts[:,2])
        return

if __name__ == '__main__':
    test = catalog_analysis("test_data/dm_files/snapshot_000","test_data/halo_files/outbgc2_0.list",out_path="./output/",down_sampling=10)
    test.build_WL_signal()

