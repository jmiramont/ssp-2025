from mcsm_benchs.benchmark_utils import MethodTemplate
import numpy as np
from numpy import pi as pi
import scipy.signal as sg
from src.stft_wrap.stft_wrapper import compute_istft, compute_spectrogram_and_zeros
from src.utilities.utilstf import get_round_window, get_spectrogram, reconstruct_signal_3
def hard_thresholding(signal, coeff=3.0, window=None, Nfft=None, dict_output=False):
    if Nfft is None:
        Nfft = 2*len(signal)
    
    if window is None:    
        window, a = get_round_window(Nfft)
    
    _, stft_whole = get_spectrogram(signal,window=window,Nfft=Nfft)

    stft_whole = stft_whole[0:Nfft//2+1,:]
    gamma = np.median(np.abs(np.real(stft_whole)))/0.6745
    thr = coeff*np.sqrt(2)*gamma

    mask = np.abs(stft_whole)
    mask[mask<thr] = 0
    mask[mask>=thr] = 1
    mask = mask.astype(bool)
    # mask[:] = 1
    # xr = np.real(reconstruct_signal_3(mask, stft_whole, window=window)) 
    xr = compute_istft(signal,nfft=Nfft,mask=mask)

    if dict_output:
        return {'xr': xr, 'mask': mask, 'stft':stft_whole} 
    else:
        return xr


class NewMethod(MethodTemplate):
    def __init__(self):
        self.id = 'thresholding_hard'
        self.task = 'denoising'


    def method(self, signal, *args, **kwargs):
        signal_output = hard_thresholding(signal, *args, **kwargs)          
        return signal_output