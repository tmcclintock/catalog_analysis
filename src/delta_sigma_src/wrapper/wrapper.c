#include "wrapper.h"

int interface(double*R,int NR,
	      double*xi_hm,
	      cosmology cosmo,
	      interface_parameters*params,
	      wrapper_output*outputs){
  int i;

  double*sigma_r=outputs->sigma_r;
  double*delta_sigma=outputs->delta_sigma;
  double*mis_sigma_r=outputs->mis_sigma_r;
  double*mis_delta_sigma=outputs->mis_delta_sigma;

  //Used to hold integration errors for some routines
  double*err=(double*)malloc(NR*sizeof(double));

  double Mass=params->Mass;
  double concentration=params->concentration;
  double Rmis=params->Rmis;
  double fmis=params->fmis;
  int delta=params->delta;

  int timing=params->timing;
  int miscentering=params->miscentering;
  int*flow_control=params->flow_control;

  double time=omp_get_wtime();
  calc_sigma_r(R,Mass,concentration,delta,R,xi_hm,NR,sigma_r,err,cosmo);
  if (timing){
    printf("sigma_r time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }
  calc_delta_sigma(R,Mass,concentration,delta,R,sigma_r,NR,delta_sigma,err,cosmo);
  if (timing){
    printf("delta_sigma time = %f\n",omp_get_wtime()-time);fflush(stdout);
    time=omp_get_wtime();
  }

  if (miscentering){
    calc_miscentered_sigma_r(R,Mass,concentration,delta,Rmis,R,sigma_r,NR,mis_sigma_r,err,cosmo);
    if (timing){
      printf("mis_sigma_r time = %f\n",omp_get_wtime()-time);fflush(stdout);
      time=omp_get_wtime();
    }
    calc_miscentered_delta_sigma(R,Mass,concentration,delta,Rmis,R,sigma_r,mis_sigma_r,NR,mis_delta_sigma,err,cosmo);
    if (timing){
      printf("mis_delta_sigma time = %f\n",omp_get_wtime()-time);fflush(stdout);
      time=omp_get_wtime();
    }
  }

  free(err);
  return 0;
}

int python_interface(double*R,int NR,double*xi_hm,
		     double h,double om,double ode,double ok,
		     double Mass,double concentration,
		     double Rmis,double fmis,
		     int delta,
		     int*flow_control,int timing,int miscentering,
		     double*sigma_r,double*delta_sigma,
		     double*mis_sigma_r,double*mis_delta_sigma){
  
  cosmology*cosmo = (cosmology*)malloc(sizeof(cosmology));
  cosmo->h = h;
  cosmo->om = om;
  cosmo->ode = ode;
  cosmo->ok = ok;

  interface_parameters*params=
    (interface_parameters*)malloc(sizeof(interface_parameters));
  params->Mass=Mass;
  params->concentration=concentration;
  params->Rmis=Rmis;
  params->fmis=fmis;
  params->delta=delta;
  params->timing=timing; //1 is true
  params->miscentering=miscentering; //1 is true

  wrapper_output*outputs=(wrapper_output*)malloc(sizeof(wrapper_output));
  outputs->sigma_r=sigma_r;
  outputs->delta_sigma=delta_sigma;
  outputs->mis_sigma_r=mis_sigma_r;
  outputs->mis_delta_sigma=mis_delta_sigma;

  interface(R,NR,xi_hm,*cosmo,params,outputs);

  return 0;
}
		     
