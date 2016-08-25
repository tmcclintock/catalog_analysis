"""
This is the python wrapper that calls the python_wrapper()
c function. This interfaces through c_types so that the user
doesn't have to.
"""
import numpy as np
import ctypes
from ctypes import c_double,c_int,POINTER,cdll

def build_Delta_Sigma(R,xi_hm,cosmo_dict,input_params):
    dslib = cdll.LoadLibrary("src/Build_Delta_Sigma_Library.so")
    interface = dslib.python_interface
    interface.restype = c_int

    """
    Arguments are: 
    R,NR,xi_hm
    h,om,ode,ok,
    Mass,concentration,
    Rmis,fmis,
    delta,
    flow_control,timing,miscentering
    sigma_r,delta_sigma,
    """

    interface.argtypes=[POINTER(c_double),c_int,\
                            POINTER(c_double),\
                            c_double,c_double,c_double,c_double,\
                            c_double,c_double,\
                            c_double,c_double,\
                            c_int,\
                            POINTER(c_int),c_int,c_int,\
                            POINTER(c_double),POINTER(c_double),\
                            POINTER(c_double),POINTER(c_double)]
    NR = len(R)
    if not (NR==len(xi_hm)):
        print "Error: len(R) != len(xi_hm)"
        return
    
    R_in = R.ctypes.data_as(POINTER(c_double))
    xi_hm_in = xi_hm.ctypes.data_as(POINTER(c_double))

    Mass,concentration,delta,timing = input_params["Mass"],input_params["concentration"],input_params["delta"],input_params["timing"]
    h,om,ode,ok = cosmo_dict['h'],cosmo_dict['om'],cosmo_dict['ode'],cosmo_dict['ok']
    miscentering = input_params["miscentering"]
    if miscentering == 1:
        Rmis,fmis = input_params["Rmis"],input_params["fmis"]
    else:
        Rmis = fmis = 0
    flow_control = np.zeros(1).ctypes.data_as(POINTER(c_int))

    sigma_r = np.zeros(NR)
    sigma_r_in = sigma_r.ctypes.data_as(POINTER(c_double))
    delta_sigma = np.zeros(NR)
    delta_sigma_in = delta_sigma.ctypes.data_as(POINTER(c_double))
    mis_sigma_r = np.zeros(NR)
    mis_sigma_r_in = mis_sigma_r.ctypes.data_as(POINTER(c_double))
    mis_delta_sigma = np.zeros(NR)
    mis_delta_sigma_in = mis_delta_sigma.ctypes.data_as(POINTER(c_double))

    result = interface(R_in,NR,xi_hm_in,\
                           h,om,ode,ok,\
                           Mass,concentration,\
                           Rmis,fmis,\
                           delta,\
                           flow_control,timing,miscentering,\
                           sigma_r_in,delta_sigma_in,\
                           mis_sigma_r_in,mis_delta_sigma_in)
    
    #Now build a dictionary and return it
    if miscentering == 0:
        return_dict = {"R":R,"xi_hm":xi_hm,"sigma_r":sigma_r,"delta_sigma":delta_sigma}
    else:
        return_dict = {"R":R,"xi_hm":xi_hm,"sigma_r":sigma_r,"delta_sigma":delta_sigma,"miscentered_sigma_r":mis_sigma_r,"miscentered_delta_sigma":mis_delta_sigma,"full_delta_sigma":(1.-fmis)*delta_sigma+fmis*mis_delta_sigma}

    return return_dict
