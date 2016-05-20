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
    TODO: also allow for passing of a TreeCorr dictionary, some kind of 
    data structure that will handle flow control, output paths for things,
    a set of mass bins, and a number of jackknife regions.
    """
    def __init__(self,catalog_name,random_name,catalog_name2=None,random_name2=None,\
                 treecorr_dict=None,flow_control=None,temp_path=None,):
        self.catalog_name  = catalog_name
        self.random_name   = random_name
        self.catalog_name2 = catalog_name2
        self.random_name2  = random_name2
        self.treecorr_dict = treecorr_dict
        self.flow_control  = flow_control
        self.temp_path     = temp_path
        return

    def analyze(self):
        """
        Parts of this function:
        1) Do exception handling. Check for file existence,
        check that the files have data in them.
        2) Interface with the functions that call TreeCorr.
        3) Evaluate flow control and figure out if 
        the WL and/or jackknifing should be calculated.
        """
        #Check that the files exist
        if not os.path.isfile(self.catalog_name):
            raise Exception("Error: catalog file not found.")
        if not os.path.isfile(self.random_name):
            raise Exception("Error: random file not found.")
        if self.catalog_name2 is not None and not os.path.isfile(self.catalog_name2):
            raise Exception("Error: catalog2 file not found.")
        if self.random_name2 is not None and not os.path.isfile(self.random_name2):
            raise Exception("Error: random2 file not found.")

        #Check that the files are non-zero
        if not os.stat(self.catalog_name).st_size > 0:
            raise Exception("Error: catalog file empty.")
        if not os.stat(self.random_name).st_size > 0:
            raise Exception("Error: random file empty.")
        if self.catalog_name2 is not None and not os.stat(self.catalog_name2).st_size > 0:
            raise Exception("Error: catalog2 file empty.")
        if self.random_name2 is not None and not os.stat(self.random_name2).st_size > 0:
            raise Exception("Error: random2 file empty.")

        #Interface with TreeCorr
        #TODO
        return

if __name__ == '__main__':
    test = catalog_analysis("test_data/test_cat.txt","test_data/test_rand.txt")
    test.analyze()
