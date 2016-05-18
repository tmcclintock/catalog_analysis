# catalog_analysis
This is an interface to a suite of applications that take in cluster, dark matter, and random catalogs and generates a variety of data vectors associated with them.

The algorithm will look like the following:

1) User supplies a cluster catalog name, a DM particle catalog name, and random particle catalog names. Optional inputs include: path to a treecorr dictionary, options to do certain operations, a specific path to put outputs and temporary files, a set of mass bins (or else a default set of bins will be used), number of jackknife regions if jackknifing is turned on.

2) Check that all files exist and are more than 0 lines long.

3) Scan the cluster file and determine minimum possible cluster that we can reliably work with, for instance one with >=100 particles.

4) Create smaller cluster catalogs broken up by mass bins.

5) Calculate requested operations including: angular correlation funciton, 3D correlation function (both from TreeCorr), and a simulated WL signal.

6) Write all current outputs to file.

7) If jackknifing was requested, then scan catalogs to determine the bounds, assuming it is square.

7.1) Create jackknife subdivision catalogs.

7.2) Calculate correlations between subdivisions and resum to get both/either the angular correlation function and 3D correlation function. Calculate the covariance matrices at this stage.

7.3) If requested, generate the WL signal for each jackknife realization and calculate the covariance matrix.
