"""
This file contains stand-alone versions of the functions
needed to interface properly with TreeCorr,
including the build_treecorr_dict()
function and the run_treecorr() function.
"""
import treecorr
import numpy as np
import pygadgetreader as pgr

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

def run_treecorr(config,dm_base_file_path,halo_base_file_path,
                 dm_random_path,halo_random_path,
                 redshift,mass_bounds,N_jks,ndivs,side_length,DS,mapping,
                 cf_output_base):
    print "Performing the treecorr analysis."
    lMmin, lMmax = mass_bounds
    step = side_length/ndivs
    
    #Read in the randoms.
    rand_dm = np.loadtxt(dm_random_path)
    rand_halo = np.loadtxt(halo_random_path)
    print rand_dm.shape, rand_halo.shape

    #Loop over all JK pairs
    for i in xrange(0,1):#ndivs):
        for j in xrange(0,1):#ndivs):
            for k in xrange(0,1):#ndivs:
                #Move the DM randoms
                rand_dm[:,0] += i*step #x
                rand_dm[:,1] += j*step #y
                rand_dm[:,2] += k*step #z

                unmapped_dm_jkindex = k + ndivs*j + ndivs*ndivs*i
                dm_jkindex = mapping[unmapped_dm_jkindex]
                dm_parts_all = pgr.readsnap(dm_base_file_path%dm_jkindex,"pos","dm",single=True,suppress=True)
                N_dm = len(dm_parts_all)
                f_keep = 1./DS
                dm_parts = []
                for dmp in range(len(dm_parts_all)):
                    r = np.random.random()
                    if r <f_keep:
                        dm_parts.append(dm_parts_all[dmp])
                dm_parts = np.array(dm_parts)
                print dm_parts_all.shape,dm_parts.shape

                #Loop over the halo indices
                for ii in xrange(0,1):#ndivs):
                    for jj in xrange(0,1):#ndivs):
                        for kk in xrange(0,1):#ndivs):
                            halo_jkindex = kk + ndivs*jj + ndivs*ndivs*ii
                            
                            #Move the halo randoms
                            rand_halo[:,0] += ii*step #x
                            rand_halo[:,1] += jj*step #y
                            rand_halo[:,2] += kk*step #z

                            halos = np.atleast_2d(np.loadtxt(halo_base_file_path%(redshift,lMmin,lMmax,redshift,lMmin,lMmax,halo_jkindex)))
                            halo_pos = halos[:,1:] #Just the positions
                            print halo_pos.shape

                            #Perform the actual treecorr analysis
                            dm_cat = treecorr.Catalog(x=dm_parts[:,0].copy(),y=dm_parts[:,1].copy(),z=dm_parts[:,2].copy(),config=config)
                            halo_cat = treecorr.Catalog(x=halo_pos[:,0].copy(),y=halo_pos[:,1].copy(),z=halo_pos[:,2].copy(),config=config)
                            random_dm_cat = treecorr.Catalog(x=rand_dm[:,0].copy(),y=rand_dm[:,1].copy(),z=rand_dm[:,2].copy(),config=config)
                            random_halo_cat = treecorr.Catalog(x=rand_halo[:,0].copy(),y=rand_halo[:,1].copy(),z=rand_halo[:,2].copy(),config=config)
                            print "\tPerforming cross correlation between DM%d and halos%d"%(dm_jkindex,halo_jkindex)
                            DD = treecorr.NNCorrelation(config)
                            DR = treecorr.NNCorrelation(config)
                            RD = treecorr.NNCorrelation(config)
                            RR = treecorr.NNCorrelation(config)
                            DD.process(dm_cat,halo_cat)
                            print "\tDD complete."
                            RR.process(random_dm_cat,random_halo_cat)
                            print "\tRR complete."
                            DR.process(dm_cat,random_dm_cat)
                            print "\tDR complete."
                            RD.process(halo_cat,random_halo_cat)
                            print "\tRD complete."
                            DD.write(cf_output_base%(redshift,lMmin,lMmax,redshift,lMmin,lMmax,dm_jkindex,halo_jkindex),RR,DR,RD)
                            print "\tWriting complete."

                            #Untranslate the halo randoms
                            rand_halo[:,0] -= ii*step #x
                            rand_halo[:,1] -= jj*step #y
                            rand_halo[:,2] -= kk*step #z
                            continue # end kk
                        continue # end jj
                    continue # end ii
                
                #Untranslate the dm randoms
                rand_dm[:,0] -= i*step #x
                rand_dm[:,1] -= j*step #y
                rand_dm[:,2] -= k*step #z
                continue # end k
            continue # end j
        continue # end i
    print "\tCompleted the treecorr analysis."
    return
