"""
This is a stand-alone file that sorts out halos
that are only within the mass bounds of interest.
"""
import numpy as np

def sort_halos(filtered_halo_path,redshift,bounded_halo_path,mass_bounds):
    print "Sorting filtered halos."
    halos = np.genfromtxt(filtered_halo_path)
    Mass,x,y,z = halos.T
    lM = np.log10(Mass)
    if mass_bounds is None:
        Mmin,Mmax = min(Mass),max(Mass)
        lMmin,lMmax = np.log10(Mmin),np.log10(Mmax)
        mass_bounds = [Mmin,Mmax]
        lmb = [lMmin,lMmax]
        print "\tCreatingmass bins between %.1f and %.1f"%(lMmin,lMmax)
    Mmin,Mmax = mass_bounds
    lMmin,lMmax = np.log10(mass_bounds)
    outfile = open(bounded_halo_path%(redshift,Mmin,Mmax),"w")
    for i in range(len(halos)):
        if Mass[i] > Mmin and Mass[i] < Mmax or Mass[i] == Mmin or Mass[i] == Mmax:
            outfile.write("%e %e %e %e\n"%(Mass[i],x[i],y[i],z[i]))
    outfile.close()
    print "\tHalos sorted into bins."
    return [mass_bounds, lmb]
