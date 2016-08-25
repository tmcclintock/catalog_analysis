/*
  This is an interface between the DeltaSigma code and 
  a piece of code that can pass in a cosmology and 
  a power spectrum.
  In addition, the next thing to be implemented will be 
  a boolean array that provides flow control with regards
  to which parts of the code to run.
*/
#include "../delta_sigma/delta_sigma.h"
#include "../sigma_r/sigma_r.h"
#include "../cosmology/cosmology.h"

#ifndef INTERFACE
#define INTERFACE
typedef struct interface_parameters{
  double Mass, concentration;
  double Rmis, fmis;
  int delta;
  int*flow_control;
  int timing, miscentering;
}interface_parameters;
#endif

#ifndef WRAPPER_OUTPUT
#define WRAPPER_OUTPUT
typedef struct wrapper_output{
  double*sigma_r;
  double*delta_sigma;
  double*mis_sigma_r;
  double*mis_delta_sigma;
}wrapper_output;
#endif

int interface(double*R,int NR,
	      double*xi_hm,
	      cosmology cosmo,
	      interface_parameters*params,
	      wrapper_output*outputs);

int python_interface(double*R,int NR,
		     double*xi_hm,
		     double h,double om,double ode,double ok,
		     double Mass,double concentration,
		     double Rmis,double fmis,
		     int delta,
		     int*flow_control,int timing,
		     int miscentering,
		     double*sigma_r,double*delta_sigma,
		     double*mis_sigma_r,double*mis_delta_sigma);

