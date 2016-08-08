"""
This is the top level object that contains the interface to the 
catalog analysis code. This code recieves a catalog and a variety 
of possible options and performs certain calculations on them
according to the user's specification.

For more information visit:
https://github.com/tmcclintock/catalog_analysis

Step 1:
The user must supply a catalog name and a random catalog name corresponding to that 
catalog. They could also provide a second catalog and random catalog if not than
a TreeCorr autocorrelation is performed and that's it. If yes then the full
WL analysis is run.
"""

import os

class catalog_analysis(object):
    """
    The object for taking in the catalog and building various quantities from it,
    including the mass function N(M), the 3D correlation function xi(r),
    and the weak lensing signal DeltaSigma.

    :dm_cat:
        The path to the dark matter (DM) catalog.

    :dm_rand:
        The path to the random catalog for the DM catalog.

    :halo_cat:
        The path to the halo catalog.

    :halo_rand:
        The path to the random catalog for the halos.

    :treecorr_dict:
        A treecorr dictionary for running TreeCorr. This contains
        all of the information need to construct xi(r).

    :outdir_path:
        A path to the top-level directory that will contain all outputs.
        This should have enough storage to store at least 20% of the DM particles
        used in the original snapshot.

    :debug:
        An option for printing out diagnostic information.
    """
    def __init__(self,dm_cat,dm_rand,halo_cat,halo_rand,\
                 treecorr_dict=None,outdir_path="./",debug=False):
        self.dm_cat  = dm_cat
        self.dm_rand   = dm_rand
        self.halo_cat = halo_cat
        self.halo_rand  = halo_rand
        self.treecorr_dict = treecorr_dict
        self.outdir_path = outdir_path
        self.debug = debug
        return

    def analyze(self):
        """
        Parts of this function:
        1. Do exception handling. Check for file existence,
        check that the files have data in them.
        2. Interface with the functions that call TreeCorr.
        """
        #Check that the files exist
        if not os.path.isfile(self.dm_cat):
            raise Exception("Error: DM catalog not found.")
        if not os.path.isfile(self.dm_rand):
            raise Exception("Error: DM random catalog not found.")
        if self.halo_cat is not None and not os.path.isfile(self.halo_cat):
            raise Exception("Error: halo catalog not found.")
        if self.halo_rand is not None and not os.path.isfile(self.halo_rand):
            raise Exception("Error: halo random catalog not found.")

        #Check that the files are non-zero
        if not os.stat(self.dm_cat).st_size > 0:
            raise Exception("Error: DM catalog file empty.")
        if not os.stat(self.dm_rand).st_size > 0:
            raise Exception("Error: DM random catalog empty.")
        if self.halo_cat is not None and not os.stat(self.halo_cat).st_size > 0:
            raise Exception("Error: halo catalog empty.")
        if self.halo_rand is not None and not os.stat(self.halo_rand).st_size > 0:
            raise Exception("Error: halo random catalog empty.")

        #TODO: Call the other objects that do the necessary parts in the algorithm

        return 0

if __name__ == '__main__':
    test = catalog_analysis("test_data/test_cat.txt","test_data/test_rand.txt",\
                                "test_data/test_cat.txt","test_data/test_rand.txt")
    if not test.analyze():
        print "Success."
    else:
        print "Failure."
