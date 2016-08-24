"""
This is a stand-alone version of the jackknife resum algorithm that takes
all of the correlation function jackknife files and resums them
into leave-one-out versions. These are the actual jackknife resampled
quantities.
"""
import numpy as np
import os,sys

def resum_cf(cf_singles_input,N_jks,ndivs,redshift,mass_bounds,cf_resum_output_base,cf_full_output_base):
    print "Resumming treecorr outputs."
    lMmin,lMmax = mass_bounds
    header = "R_nom\tmeanR\tmeanlogR\txi\tsigmaxi\tDD\tRR\tDR\tRD\tnpairs\n"

    #Read in a single input file to get the data shape
    data = np.genfromtxt(cf_singles_input%(redshift,lMmin,lMmax,redshift,lMmin,lMmax,0,0))
    R_full,R_mean_full,ln_mean_R_full,xi_full,sigmaxi_full,DD_full,RR_full,DR_full,RD_full,npairs_full = data.T
    xi_full = np.fabs(0*xi_full)
    DD_full = np.fabs(0*DD_full)
    RR_full = np.fabs(0*RR_full)
    DR_full = np.fabs(0*DR_full)
    RD_full = np.fabs(0*RD_full)
    
    #Loop over the jackknife_indices to get the full versions
    for i in range(0,1):#N_jks):
        for j in range(0,1):#N_jks):
            in_data = np.genfromtxt(cf_singles_input%(redshift,lMmin,lMmax,redshift,lMmin,lMmax,i,j))
            R,R_mean,ln_mean_R,xi,sigmaxi,DD,RR,DR,RD,npairs = in_data.T
            DD_full += DD
            RR_full += RR
            DR_full += DR
            RD_full += RD
    xi_full = (DD_full-DR_full-RD_full+RR_full)/RR_full
    out_full = np.array([R_full,DD_full,RR_full,DR_full,RD_full,xi_full]).T
    np.savetxt(cf_full_output_base%(redshift,lMmin,lMmax,redshift,lMmin,lMmax),out_full)

    #Loop over the jackknife indices to create the leave-one-out resamplings
    for jkindex in range(0,1):#N_jks):
        xi_LOO = np.fabs(xi_full*0)
        DD_LOO = DD_full.copy()
        RR_LOO = RR_full.copy()
        DR_LOO = DR_full.copy()
        RD_LOO = RD_full.copy()
        #First leave out the DM particle indices
        for j in range(0,1):#N_jks):
            LOO_data = np.genfromtxt(cf_singles_input%(redshift,lMmin,lMmax,redshift,lMmin,lMmax,j,jkindex))
            R,R_mean,ln_mean_R,xi,sigmaxi,DD,RR,DR,RD,npairs = LOO_data.T
            DD_LOO -= DD
            RR_LOO -= RR
            DR_LOO -= DR
            RD_LOO -= RD
            continue # end j
        #Now leave out the halo particle indices
        #Note: jkindex,jkindex has already been left out
        for j in range(0,1):#N_jks):
            if j == jkindex: continue
            LOO_data = np.genfromtxt(cf_singles_input%(redshift,lMmin,lMmax,redshift,lMmin,lMmax,jkindex,j))
            R,R_mean,ln_mean_R,xi,sigmaxi,DD,RR,DR,RD,npairs = LOO_data.T
            DD_LOO -= DD
            RR_LOO -= RR
            DR_LOO -= DR
            RD_LOO -= RD
            continue # end j
        #The next line throws an error during testing since we don't have
        #more than one jackknife region.
        xi_LOO = (DD_LOO-DR_LOO-RD_LOO+RR_LOO)/(RR_LOO)
        out_LOO = np.array([R_full,DD_LOO,RR_LOO,DR_LOO,RD_LOO,xi_LOO]).T
        np.savetxt(cf_resum_output_base%(redshift,lMmin,lMmax,redshift,lMmin,lMmax,jkindex),out_LOO)
        continue # end jkindex

    import matplotlib.pyplot as plt
    plt.loglog(R_full,xi_full)
    plt.show()

    print "\tResummation complete."
    return
