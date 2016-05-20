#! /usr/bin/env python

# Correlation function calculation
## User supplies data and random RA and DEC (and optionally R) and TreeCorr config dict
## Function returns correlation function (xi) and logr
def calculate_xi(data_ra, data_dec, rand_ra, rand_dec, config, data_r=None, rand_r=None) :
	## Imports
	import numpy as np
	import treecorr
	
	## Set up TreeCorr catalog objects
	if data_r is not None and rand_r is not None :
		dcat = treecorr.Catalog(ra=data_ra, dec=data_dec, r=data_r, ra_units='deg', dec_units='deg')
		rcat = treecorr.Catalog(ra=rand_ra, dec=rand_dec, r=rand_r, ra_units='deg', dec_units='deg')
	else :
		dcat = treecorr.Catalog(ra=data_ra, dec=data_dec, ra_units='deg', dec_units='deg')
		rcat = treecorr.Catalog(ra=rand_ra, dec=rand_dec, ra_units='deg', dec_units='deg')
	
	## Run TreeCorr processes (i.e. do pair counts)
	### Data-Data
	dd = treecorr.NNCorrelation(config)
	dd.process(dcat)
	### Data-Random
	dr = treecorr.NNCorrelation(config)
	dr.process(dcat, rcat)
	### Random-Random
	rr = treecorr.NNCorrelation(config)
	rr.process(rcat)
	
	## Calculate CF
	xi = dd.calculateXi(rr, dr)[0]
	
	## Return values
	return (xi, dd.logr)

# Great Circle distance converter
## Converts chord length to true great circle physical or angular distances
## User provides chord lengths C in radians (and optionally line of sight distance R)
## Function returns great circle distance (in same units as R or in radians)
def great_circle_dist(C, R=1) :
	## Imports
	import numpy as np
	
	## Calculate
	d = 2.*R*np.arcsin(C/2.)
	
	## Return distances
	return (d)