"""
This is a stand-alone file that contains a method for 
re-ordering the jackknife files to find the correct mapping
from the LGADGET jackknifing numbering system
to our numbering system.
"""
import numpy as np
import pygadgetreader as pgr

def reorder_jackknifes(N_jk,side_length,dm_path):
    print "Finding the mapping of LGADGET jackknifes to our jackknifes."
    ndivs = int(round(N_jk**(1./3.)))
    step = side_length/ndivs
    mapping = {}
    for i in xrange(0,2): #N_jk):
        pos = pgr.readsnap(dm_path+".%d"%i,"pos","dm",single=True,suppress=True)
        xi,yi,zi = [int(q) for q in np.mean(pos,0)/step]
        jkindex = zi + ndivs*yi + ndivs*ndivs*xi
        mapping[jkindex] = i
    print "\tMapping complete on %d files."%N_jk
    return mapping


"""
def reorder_jackknifes(self):
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
"""
