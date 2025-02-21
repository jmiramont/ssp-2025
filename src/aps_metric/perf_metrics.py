# One possible metric of musical noise:
import librosa
import numpy as np
import scipy.signal as sg
from numpy.linalg import norm
import os

# Define a performance function from Matlab
from mcsm_benchs.MatlabInterface import MatlabInterface

def aps_measure(x,n,x_hat,fs=None,*args,**kwargs):
    # APS from PEASS matlab implementation
    paths = [os.path.join('src','aps_metric','peass'),
            os.path.join('..','src','aps_metric','peass')
            ]
    
    mlint = MatlabInterface('APS_wrapper', add2path=paths) 
    aps_wrapper = mlint.matlab_function # A python function handler to the method.
    if fs is None:
        # fs = 11025
        fs = 8000
    
    aps = aps_wrapper(x,n,x_hat,fs)
    
    return aps

def musical_noise_measure_aps(x,x_hat,fs=None,**kwargs):
    if 'scaled_noise' in kwargs.keys():
        n = kwargs['scaled_noise']

    APS = aps_measure(x,n,x_hat,fs=fs) 

    return APS

from oct2py import octave
def aps_measure_octave(x,n,x_hat,fs=None,*args,**kwargs):

    if fs is None:
        # fs = 11025
        fs = 8000

    # APS from PEASS matlab implementation
    paths = [os.path.join('src','aps_metric','peass'),
            os.path.join('src','aps_metric','peass','gammatone'),
            os.path.join('..','src','aps_metric','peass'),
            os.path.join('..','src','aps_metric','peass','gammatone'),
            ]
    _ = octave.addpath(paths[0])
    _ = octave.addpath(paths[1])

    aps = octave.APS_wrapper_8(x,n,x_hat,fs,nout=1)
    return aps

def musical_noise_measure_aps_8(x,x_hat,fs=None,**kwargs):
    if 'scaled_noise' in kwargs.keys():
        n = kwargs['scaled_noise']

    APS = aps_measure_octave(x,n,x_hat,fs=fs) 

    return APS