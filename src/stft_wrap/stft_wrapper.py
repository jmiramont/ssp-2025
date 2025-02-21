import numpy as np
import ctypes
import os
import matplotlib.pyplot as plt

# Load the shared library
if os.name == "nt":
    stft_lib = ctypes.CDLL("./stft.dll")
else:
    stft_lib = ctypes.CDLL("./src/stft_wrap/libstft.so")

stft_lib.compute_spectrogram_and_zeros.argtypes = [
    ctypes.POINTER(ctypes.c_double), 
    ctypes.c_int, 
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]

stft_lib.compute_istft.argtypes = [
    ctypes.POINTER(ctypes.c_double), 
    ctypes.c_int,
    ctypes.c_int, 
    ctypes.POINTER(ctypes.c_bool), 
    ctypes.POINTER(ctypes.c_double)
]

# Python wrapper function for STFT
def compute_spectrogram_and_zeros(signal):

    signal = np.ascontiguousarray(signal, dtype=np.float64)
    original_length = len(signal)
    nfft = 2*original_length

    signal_aux = np.zeros(nfft,dtype=np.float64)
    signal_aux[original_length//2:nfft-original_length//2]=signal 

    signal_aux = np.ascontiguousarray(signal_aux)
    signal_length = len(signal_aux)

    # spectrogram_out = np.zeros((original_length,nfft//2+1), dtype=np.float64)
    zeros_x = np.zeros((nfft,), dtype=np.int32)
    zeros_y = np.zeros((nfft,), dtype=np.int32)
    zeros_count = np.zeros((1,), dtype=np.int32)

    stft_lib.compute_spectrogram_and_zeros(
        signal_aux.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        np.int32(signal_length), 
        np.int32(nfft),
        zeros_x.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        zeros_y.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        zeros_count.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )

    zeros = np.zeros((zeros_count[0],2))
    zeros[:,1] = zeros_x[0:zeros_count[0]]
    zeros[:,0] = zeros_y[0:zeros_count[0]]
    # return real_out[:num_frames.value, :], imag_out[:num_frames.value, :]
    # spectrogram_out[original_length//2-N:int(3*original_length//2)-N,:]
    return zeros

# Python wrapper function for ISTFT
def compute_istft(signal, nfft=None, mask=None):
    signal = np.ascontiguousarray(signal, dtype=np.float64)
    original_length = len(signal)
    if nfft is None:
        nfft = 2*original_length

    prec = 1e-15
    L = np.sqrt(nfft)
    l = np.floor(np.sqrt(-nfft*np.log(prec)/np.pi))+1
    N = int(l+1)

    signal_aux = np.zeros(2*original_length, dtype=np.float64)
    signal_aux[original_length//2:original_length//2+original_length]=signal 

    signal_length = len(signal_aux)


    if mask is None:
        mask = np.ones((original_length,nfft),dtype=bool)


    reconstructed_signal = np.zeros(signal_length, dtype=np.float64)

    stft_lib.compute_istft(
        signal_aux.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        signal_length,
        nfft,
        mask.T.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
        reconstructed_signal.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    return reconstructed_signal [N:N+original_length]

# Test the wrapper
if __name__ == "__main__":
    fs = 2**13 # 16 kHz
    duration = 1.0
    t = np.linspace(0, duration, int(fs * duration), endpoint=False)
    # signal = np.sin(2 * np.pi * 32 * t)
    # signal = np.random.randn(fs) +  30*np.sin(2 * np.pi * 128 * t)

    import librosa
    signal, fs = librosa.load(librosa.ex('trumpet'), duration=1.5, sr=8000)
    signal = signal[0:fs]
    # signal[128] = 256

    # real_part, imag_part, window = compute_stft(signal)
    # zeros = compute_spectrogram_and_zeros(signal)

    # print(spectrogram.shape)
    r_signal = compute_istft(signal, nfft=2048)

    plt.plot(signal)
    plt.plot(r_signal,'--r', alpha=0.5)

    # plt.plot(zeros[:,1],zeros[:,0],'r.')
    # plt.axis('square')
    # plt.colorbar()

    # plt.plot(window)    
    
    plt.show()

    # reconstructed_signal = compute_istft(real_part, imag_part)

    # print("Original Signal First 10 Samples:", signal[5000:5010])
    # print("Reconstructed Signal First 10 Samples:", reconstructed_signal[5000:5010])

