"""
This is a stand-alone version of the jackknife_halos function
that was originally within the catalog_analysis.py file.
"""
import numpy as np
import os

def jackknife_halos(side_length,ndivs,redshift,log10_mass_bounds,
                    bounded_halo_path,
                    bounded_jk_output_directory_base,
                    jk_halo_filename):
    """
    side_length: length of the simulation side

    ndivs: number of divisions along a side. Usually 2 or 4 or 8

    redshift: redshift of the snapshot

    log10_mass_bounds: the mass bounds of these halos, in log base 10

    bounded_halo_path: the path+filename for the bounded halos
    that isn't jackknifed yet. This is unfilled formatted.

    bounded_jk_output_directory_base: the unfilled formatted path to the 
    directory that will contain the jackknifed files for this mass bound

    jk_halo_filename: the filename for the jackknifed halo files
    for this mass bound. This is unfilled formatted.

    """
    lMmin,lMmax = log10_mass_bounds
    print "Jackknifing the halo file with mass bounds of %.1f to %.1f"%(lMmin,lMmax)
    #First create the output directory for the jackknifed files
    os.system("mkdir -p %s"%bounded_jk_output_directory_base%(redshift,lMmin,lMmax))
    step = side_length/ndivs
    halos = np.genfromtxt(bounded_halo_path%(redshift,lMmin,lMmax))
    M,x,y,z = halos.T
    xi = np.floor(x/step)
    yi = np.floor(y/step)
    zi = np.floor(z/step)
    jkindices = zi + ndivs*yi + ndivs*ndivs*xi

    outlist = []
    for jk in range(ndivs**3):
        #print bounded_jk_output_directory_base
        #print jk_halo_filename
        outlist.append(open(bounded_jk_output_directory_base%(redshift,lMmin,lMmax)+jk_halo_filename%(redshift,lMmin,lMmax,jk),"w"))
    
    for jki in range(len(jkindices)):
        outlist[int(jkindices[jki])].write("%e %e %e %e\n"%(M[jki],x[jki],y[jki],z[jki]))

    for jk in range(ndivs**3):
        outlist[jk].close()
    print "\tHalos are now jackknifed."
