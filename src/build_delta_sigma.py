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

import sys
#sys.path.insert(0,"src/delta_sigma_src/wrapper/"
#import Build_Delta_Sigma

