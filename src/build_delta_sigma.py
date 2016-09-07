"""
This is an interface file between the catalog_analysis code 
and the Build-Delta-Sigma code that is reproduced in a subdirectory located
here.

The repository for Build-Delta-Sigma can be found here:
https://github.com/tmcclintock/Build-Delta-Sigma

This code takes in the halo catalogs and the paths to the 
correlation functions. It then computes the average masses
of the halos for each jackknife-region (and the entire sim)
and run the Build-Delta-Sigma code after calling
collossus for an M-c relation.
"""

import numpy as np
import sys
sys.path.insert(0,"src/delta_sigma_src/wrapper/")
import Build_Delta_Sigma_Interface

def calculate_concentration(redshift,M200m,cosmology):
    from colossus.halo import concentration
    from colossus.cosmology import cosmology as col_cosmology
    params = {'flat': True, 'H0': cosmology['h'], 'Om0': cosmology['om'],
              'Ob0': 1.0-cosmology['om'],
              'sigma8': cosmology['sigma8'], 'ns': cosmology['ns']}
    col_cosmology.addCosmology('fiducial_cosmology',params)
    col_cosmology.setCosmology('fiducial_cosmology')
    concentration = concentration.concentration(M200m,'200m',redshift,model='diemer15')
    print "concentration = ",concentration,2.0
    return 2.0#concentration

def build_delta_sigma(cosmology,redshift,mass_bounds,halo_full_base,\
                      halo_jk_base,cf_full_base,cf_resum_base,\
                      delta_sigma_resum_output_base,\
                      delta_sigma_full_output_base):
    print "Building all delta sigma curves."
    lMmin,lMmax = mass_bounds

    #First read in the full halo catalogs
    M,Rs,x,y,z = np.genfromtxt(halo_full_base%(redshift,lMmin,lMmax)).T
    concentration = calculate_concentration(redshift,np.mean(M),cosmology)
    mean_mass_full = np.mean(M)
    mean_Rs_full = np.mean(Rs)
    mean_conc_full = np.mean(concentration)

    #Now read in the full CF file
    R_full,DD,RR,DR,RD,xi_full = np.genfromtxt(cf_full_base%(redshift,lMmin,lMmax,redshift,lMmin,lMmax)).T
    
    
    #Create a dictionary with the input parameters
    input_params = {"Mass":mean_mass_full,"delta":200,\
                    "concentration":mean_conc_full,\
                    "timing":1,"miscentering":0,}
    return_dict_full = Build_Delta_Sigma_Interface.build_Delta_Sigma(R_full.copy(),xi_full.copy(),cosmology,input_params)
    print return_dict_full.keys()

    import matplotlib.pyplot as plt
    R = return_dict_full['R']
    xi = return_dict_full['xi_hm']
    sigma_r = return_dict_full['sigma_r']
    delta_sigma = return_dict_full['delta_sigma']
    plt.loglog(R,xi)
    plt.loglog(R,sigma_r)
    plt.loglog(R,delta_sigma)
    plt.show()


