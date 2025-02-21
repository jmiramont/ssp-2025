""" This file contains a number of utilities for time-frequency analysis. 
Some functions has been modified from the supplementary code of:
Bardenet, R., Flamant, J., & Chainais, P. (2020). "On the zeros of the spectrogram of 
white noise." Applied and Computational Harmonic Analysis, 48(2), 682-705.

Those functions are:
- getSpectrogram(signal)
- findCenterEmptyBalls(Sww, pos_exp, radi_seg=1)
- getConvexHull(Sww, pos_exp, empty_mask, radi_expand=0.5)
- reconstructionSignal(hull_d, stft)
"""

import numpy as np
import scipy.signal as sg
from scipy.fft import fft, ifft
from math import factorial
from numpy import complex128, dtype, pi as pi
import pandas as pd
import os


def get_gauss_window(Nfft,L,prec=1e-6):
    l=np.floor(np.sqrt(-Nfft*np.log(prec)/pi))+1
    N = 2*l+1
    t0 = l+1
    tmt0=np.arange(0,N)-t0
    g = np.exp(-(tmt0/L)**2 * pi)    
    g=g/np.linalg.norm(g)
    return g

def get_round_window(Nfft, prec = 1e-6):
    """ Generates a round Gaussian window, i.e. same essential support in time and 
    frequency: g(n) = exp(-pi*(n/T)^2) for computing the Short-Time Fourier Transform.
    
    Args:
        Nfft: Number of samples of the desired fft.

    Returns:
        g (ndarray): A round Gaussian window.
        T (float): The scale of the Gaussian window (T = sqrt(Nfft))
    """
    # analysis window
    L=np.sqrt(Nfft)
    l=np.floor(np.sqrt(-Nfft*np.log(prec)/pi))+1
 
    N = 2*l+1
    t0 = l+1
    tmt0=np.arange(0,N)-t0
    g = np.exp(-(tmt0/L)**2 * pi)    
    g=g/np.linalg.norm(g)
    return g, L

def get_stft(x,window=None,t=None,Nfft=None):
    xrow = len(x)

    if t is None:
        t = np.arange(0,xrow)

    if Nfft is None:
        Nfft = 2*xrow

    if window is None:
        window = get_round_window(Nfft)

    tcol = len(t)
    hlength=np.floor(Nfft/4)
    hlength=int(hlength+1-np.remainder(hlength,2))

    hrow=len(window)
    
    assert np.remainder(hrow,2) == 1

    Lh=(hrow-1)//2
    tfr= np.zeros((Nfft,tcol))   
    for icol in range(0,tcol):
        ti= t[icol]; 
        tau=np.arange(-np.min([np.round(Nfft/2),Lh,ti]),np.min([np.round(Nfft/2),Lh,xrow-ti])).astype(int)
        indices= np.remainder(Nfft+tau,Nfft).astype(int); 
        tfr[indices,icol]=x[ti+tau]*np.conj(window[Lh+tau])/np.linalg.norm(window[Lh+tau])
    
    tfr=fft(tfr, axis=0) 
    return tfr

def get_istft(tfr,window=None,t=None):
    
    N,NbPoints = tfr.shape
    tcol = len(t)
    hrow = len(window) 
    Lh=(hrow-1)//2
    window=window/np.linalg.norm(window)
    tfr=ifft(tfr,axis=0)

    x=np.zeros((tcol,),dtype=complex)
    for icol in range(0,tcol):
        valuestj=np.arange(np.max([1,icol-N/2,icol-Lh]),np.min([tcol,icol+N/2,icol+Lh])).astype(int)
        for tj in valuestj:
            tau=icol-tj 
            indices= np.remainder(N+tau,N).astype(int)
            x[icol]=x[icol]+tfr[indices,tj]*window[Lh+tau]
        
        x[icol]=x[icol]/np.sum(np.abs(window[Lh+icol-valuestj])**2)
    return x


# def get_stft(signal, window = None, overlap = None):
#     """ Compute the STFT of the signal. Signal is padded with zeros.
#     The outputs corresponds to the STFT with the regular size and also the
#     zero padded version. The signal is zero padded to alleviate border effects.

#     Args:
#         signal (ndarray): The signal to analyse.
#         window (ndarray, optional): The window to use. If None, uses a rounded Gaussian
#         window. Defaults to None.

#     Returns:
#         stft(ndarray): Returns de stft of the signal.
#         stft_padded(ndarray): Returns the stft of the zero-padded signal.
#         Npad(int): Number of zeros padded on each side of the signal.
#     """
    
#     N = np.max(signal.shape)
#     if window is None:
#         window, _ = get_round_window(N)

#     Npad = N//2 #0 #int(np.sqrt(N))
#     Nfft = len(window)
    
#     if overlap is None:
#         overlap = Nfft-1
    
#     if signal.dtype == complex128:
#         signal_pad = np.zeros(N+2*Npad, dtype=complex128)
#     else:
#         signal_pad = np.zeros(N+2*Npad)

#     # signal_pad = np.zeros(N+2*Npad)
#     signal_pad[Npad:Npad+N] = signal

#     # computing STFT
#     _, _, stft_padded = sg.stft(signal_pad, window=window, nperseg=Nfft, noverlap = overlap)
    
#     # if signal.dtype == complex128:
#     #     stft_padded = stft_padded[0:Nfft//2+1,:]
        
#     stft = stft_padded[:,Npad:Npad+N]
#     return stft, stft_padded, Npad


def get_spectrogram(signal,window=None,Nfft=None,t=None,onesided=True):
    """
    Get the round spectrogram of the signal computed with a given window. 
    
    Args:
        signal(ndarray): A vector with the signal to analyse.

    Returns:
        S(ndarray): Spectrogram of the signal.
        stft: Short-time Fourier transform of the signal.
        stft_padded: Short-time Fourier transform of the padded signal.
        Npad: Number of zeros added in the zero-padding process.
    """

    N = np.max(signal.shape)
    if Nfft is None:
        Nfft = 2*N

    if window is None:
        window, _ = get_round_window(Nfft)

    if t is None:
        t = np.arange(0,N)
        
    stft=get_stft(signal,window=window,t=t,Nfft=Nfft)

    if onesided:
        S = np.abs(stft[0:Nfft//2+1,:])**2
    else:
        S = np.abs(stft)**2                
    return S, stft

def reconstruct_signal_2(mask, stft, Npad, Nfft=None, window=None, overlap=None):
    """Reconstruction using a mask given as parameter

    Args:
        mask (_type_): _description_
        stft (_type_): _description_
        Npad (_type_): _description_
        Nfft (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """

    Ni = mask.shape[1]
    if Nfft is None:
        Nfft = Ni
    
    if window is None:
        window = sg.gaussian(Nfft, np.sqrt((Nfft)/2/np.pi))
        window = window/window.sum()
    
    Nfft = len(window)
    
    if overlap is None:
        overlap=Nfft-1


    # reconstruction        
    mask_aux = np.zeros(stft.shape)
    mask_aux[:,Npad:Npad+Ni] = mask
    # t, xorigin = sg.istft(stft, window=g,  nperseg=Nfft, noverlap=Nfft-1)
    t, xr = sg.istft(mask_aux*stft, window=window, nperseg=Nfft, noverlap = overlap)
    xr = xr[Npad:Npad+Ni]
    return xr, t

def reconstruct_signal_3(mask, stft, window=None):
    mask_complete = np.zeros_like(stft, dtype=float)
    mask_complete[0:mask.shape[0],:] = mask
    mask_complete[mask.shape[0]::,:] = mask[-2:0:-1,:]
    xr = get_istft(stft*mask_complete,window=window,t=np.arange(0,stft.shape[1]))
    return xr
    

# def voss(nrows, ncols=16):
#     """Generates pink noise using the Voss-McCartney algorithm.
    
#     nrows: number of values to generate
#     rcols: number of random sources to add
    
#     returns: NumPy array
#     """
#     array = np.empty((nrows, ncols))
#     array.fill(np.nan)
#     array[0, :] = np.random.random(ncols)
#     array[:, 0] = np.random.random(nrows)
    
#     # the total number of changes is nrows
#     n = nrows
#     cols = np.random.geometric(0.5, n)
#     cols[cols >= ncols] = 0
#     rows = np.random.randint(nrows, size=n)
#     array[rows, cols] = np.random.random(n)

#     df = pd.DataFrame(array)
#     df.fillna(method='ffill', axis=0, inplace=True)
#     total = df.sum(axis=1)

#     return total.values - np.mean(total.values)

# import librosa
# import src.utilities.utilstf as utilstf

# def babble_noise(N,):
#     """Generates pink noise using the Voss-McCartney algorithm.
    
#     nrows: number of values to generate
#     rcols: number of random sources to add
    
#     returns: NumPy array
#     """
#     # Load the noise and define the mixture
#     path = os.path.dirname(utilstf.__file__)
#     babble,fs = librosa.load(os.path.join(path,'babble.wav'), sr=8000)
#     Nb = len(babble)
#     assert N<Nb, 'Cannot generate babble noise of length {}'.format(N)

#     tinit = np.random.randint(0,Nb-N-1)

#     return babble[tinit:tinit+N] - np.mean(babble[tinit:tinit+N])