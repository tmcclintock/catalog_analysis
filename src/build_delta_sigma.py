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

G = 4.517e-48 #Newton's G in Mpc^3/s^2/Solar Mass
Mpcperkm = 3.241e-20 #Mpc/km used to convert H0 to per seconds
delta_c = 1.686 #Critical collapse density

#CHANGE THIS FUNCTION TO CALL COLOSSUS
def calculate_concentration(Rs_array,M_array,cosmology):
    #for i in range(len(Rs_array)):
    #    print Rs_array[i],M_array[i]
    delta = 200.0
    h,om = cosmology['h'],cosmology['om']
    H0 = h*100.0
    rhom = om*3.*(H0*Mpcperkm*H0*Mpcperkm)/(8.*np.pi*G)/(h*h)
    Rdelta_array = (M_array/(4./3.*np.pi*rhom*delta))**(1./3.)
    return Rdelta_array/Rs_array #concentration_array

def build_delta_sigma(cosmology,redshift,mass_bounds,halo_full_base,\
                      halo_jk_base,cf_full_base,cf_resum_base,\
                      delta_sigma_resum_output_base,\
                      delta_sigma_full_output_base):
    print "Building all delta sigma curves."
    lMmin,lMmax = mass_bounds

    #First read in the full halo catalogs
    M,Rs,x,y,z = np.genfromtxt(halo_full_base%(redshift,lMmin,lMmax)).T
    concentration = calculate_concentration(Rs,M,cosmology)
    mean_mass_full = np.mean(M)
    mean_Rs_full = np.mean(Rs)
    mean_conc_full = np.mean(concentration)

    #Now read in the full CF file
    R_full,DD,RR,DR,RD,xi_full = np.genfromtxt(cf_full_base%(redshift,lMmin,lMmax,redshift,lMmin,lMmax)).T
    
    #Create a dictionary with the input parameters
    input_params = {"Mass":mean_mass_full,"delta":200,\
                    "concentration":mean_conc_full,\
                    "timing":0,"miscentering":0,}

    return_dict_full = Build_Delta_Sigma_Interface.build_Delta_Sigma(R_full,xi_full,cosmology,input_params)
    print return_dict.keys()
