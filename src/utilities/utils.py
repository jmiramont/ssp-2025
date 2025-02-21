import numpy as np
import scipy.signal as sg
import gudhi as gd
from src.stft_wrap.stft_wrapper import compute_spectrogram_and_zeros
import os
import pickle

def find_zeros_of_spectrogram(S):
    """ Find the zeros of the spectrogram S as minima in 3x3 grids.
    Includes first/last rows/columns.
    """
    aux_S = np.zeros((S.shape[0]+2,S.shape[1]+2))+np.Inf
    aux_S[1:-1,1:-1] = S
    S = aux_S
    aux_ceros = ((S <= np.roll(S,  1, 0)) &
            (S <= np.roll(S, -1, 0)) &
            (S <= np.roll(S,  1, 1)) &
            (S <= np.roll(S, -1, 1)) &
            (S <= np.roll(S, [-1, -1], [0,1])) &
            (S <= np.roll(S, [1, 1], [0,1])) &
            (S <= np.roll(S, [-1, 1], [0,1])) &
            (S <= np.roll(S, [1, -1], [0,1])) 
            )
    [y, x] = np.where(aux_ceros==True)
    pos = np.zeros((len(x), 2)) # Position of zeros in norm. coords.
    pos[:, 0] = y-1
    pos[:, 1] = x-1
    return pos

def one_sided_spectrogram_and_zeros(signal, Nfft=None, return_stft=False):
    """
    A one sided spectrogram and its zeros (for real-valued signals only)
    """
    assert np.isreal(np.all(signal)), "The signal should be real."

    N = len(signal)
    if Nfft is None:
        Nfft = 2*N
    
    window = sg.windows.gaussian(Nfft, np.sqrt(Nfft/2/np.pi))
    # print(np.sum(window**2))
    # window = window/np.sum(window**2)**0.5

    if Nfft > N:
        sigaux = np.zeros((Nfft,),dtype=signal.dtype)
        sigaux[0:N] = signal
    else:
        sigaux = signal
        
    # Compute the STFT
    frequencies, times, stft = sg.stft(sigaux,
                                    window=window, 
                                    nperseg=Nfft, 
                                    noverlap=Nfft-1, 
                                    return_onesided=True,
                                    scaling='psd'
                                    )

    stft = stft[:,0:N]

    if return_stft:
        return stft

    spectrogram = np.abs(stft)**2
    # Find zeros
    # zeros = find_local_minima(spectrogram)
    # zeros = np.array(zeros,dtype=float).T
    zeros = find_zeros_of_spectrogram(spectrogram) # This includes the zeros in the border of the plane.

    return spectrogram, zeros, stft

def invert_one_sided_stft(stft,Nfft,mask=None,smooth=False):
    """ Invert STFT computed with the function one_sided_spectrogram_and_zeros(...) 
    A filtering mask can be provided.
    """
    
    if mask is None:
        mask = np.ones_like(stft)

    N = stft.shape[1]

    window = sg.windows.gaussian(Nfft, np.sqrt(Nfft/2/np.pi))
    window = window/window.sum()

    if Nfft > N:
        stftaux = np.zeros((stft.shape[0],Nfft), dtype=complex)
        stftaux[:,0:N] = stft*mask.astype(complex)
    else:
        stftaux = np.zeros((stft.shape[0],stft.shape[1]+1), dtype=complex)
        stftaux[:stft.shape[0],:stft.shape[1]] = stft*mask.astype(complex)
        
    # Compute the STFT
    t, xr = sg.istft(stftaux,
                window=window, 
                nperseg=Nfft, 
                noverlap=Nfft-1,
                scaling='psd'
                )

    xr = xr[0:N]
    return xr

def rotated_pd(pairs):
    rotated_pairs = []
    for pair in pairs:
        rotated_pair = [(pair[0]+pair[1])/2,  pair[1]-pair[0]]
        rotated_pairs.append(rotated_pair)
    return np.array(rotated_pairs)

def compute_apf(rotated_pairs, m=0):   
    idx = rotated_pairs[:,0]<m
    return np.sum(rotated_pairs[idx,1])

def apf_mc_test(signal,
                nsim=199, 
                alpha=0.05, 
                pnorm=2.0, 
                m_max=None, 
                Nfft=None, 
                hdim=1,
                return_all=False,
                path=None):

    N = len(signal)
    if Nfft is None:
        Nfft = 2*N


    # If the simulated ppp are present, use them, otherwise, generate the simulation and
    # save it so that it can be used later (this saves time, avoiding re-simulation).
    if path is None:
        path = '.'

    # Compute APF for signal
    if N>4095:
        zeros = compute_spectrogram_and_zeros(signal)
        S = None
    else:
        S, zeros, _ = one_sided_spectrogram_and_zeros(signal, Nfft=Nfft)
    
    points = zeros[:,[1,0]]/Nfft**0.5
   
    # ac_sig = gd.AlphaComplex(points,)
    # self._st =  # gudhi promotes storage of complexes using simplex trees
      
    ac_sig = gd.AlphaComplex(points.copy())
    diagram = ac_sig.create_simplex_tree().persistence()
    # ac_sig.compute_persistent_homology() # compute persistence diagram

    pairs = [p[1] for p in diagram if p[0]==hdim]

    rotated_pairs = rotated_pd(pairs)
    apf_signal = []
    ms = np.unique(rotated_pairs[:,0])

    if m_max is None:
        m_max = np.max(ms)

    for m in ms: # Use this ms to compute apf from noise later.
        apf_signal.append(compute_apf(rotated_pairs,m=m))

    apf_signal = np.array(apf_signal)

    try:
        with open(os.path.join(path,'ppp_APF_simulations_N_{}_Nsim_{}.mcsim'.format(N,nsim)), 'rb') as handle:
            ppp_simulation = pickle.load(handle)
        rotated_pairs = ppp_simulation['rotated_pairs']
    except:
        print('MC Simulation running...')
        # Compute APF Monte Carlo simulations
        rotated_pairs = []
        for n in range(nsim):
            # Get new noise
            noise = np.random.randn(N,)
            # Compute PD and RPD

            if N>4095:
                zeros = compute_spectrogram_and_zeros(noise)
            else:
                _, zeros, _ = one_sided_spectrogram_and_zeros(noise, Nfft=Nfft)

            ac_noise = gd.AlphaComplex(zeros[:,[1,0]]/Nfft**0.5)
            diagram_noise = ac_noise.create_simplex_tree().persistence()
            pairs = [p[1] for p in diagram_noise if p[0]==hdim]
            rotated_pairs.append(rotated_pd(pairs))

        ppp_simulation = {'rotated_pairs':rotated_pairs, 'ms':ms, 'N':N, 'nsim':nsim}

        # Save simulations for next time.
        with open(os.path.join(path,'ppp_APF_simulations_N_{}_Nsim_{}.mcsim'.format(N,nsim)), 'wb') as handle:
            pickle.dump(ppp_simulation, handle, protocol=pickle.HIGHEST_PROTOCOL) 
    
    apf_noise = np.zeros((nsim,len(ms)))
    for i in range(nsim):
        for j,m in enumerate(ms):
            apf_noise[i,j] = compute_apf(rotated_pairs[i],m=m)

    # Compute test statistics 
    apf_mean = np.vstack((apf_noise,apf_signal))
    apf_mean = np.mean(apf_mean,axis=0) 

    t0 = np.linalg.norm(apf_signal[ms<=m_max]-apf_mean[ms<=m_max],ord=pnorm)
    tj = np.zeros((nsim,))
    for j in range(nsim): 
        tj[j] = np.linalg.norm(apf_noise[j][ms<=m_max]-apf_mean[ms<=m_max],ord=pnorm)

    alpha = 0.05
    k = int(alpha*(nsim+1))
    tj_sorted = np.sort(tj)[::-1]
    
    rejectH0 = t0>tj_sorted[k]

    # p-value
    p_value = np.sum(tj>t0)/nsim
    # print(p_value, p_value<alpha)
    if return_all:
        output_dict = {'apf_obs':apf_signal,
                       'apf_noise':apf_noise,
                       'apf_mean':apf_mean,
                       'rejectH0':rejectH0,
                       'p_value':p_value,
                       'k':k,
                       'ms':ms,
                       'nsims':nsim,
                       'spectrogram':S,
                       }
        return output_dict
    return rejectH0, p_value

import librosa
import src.utilities.utils as utilstf

def babble_noise(N,):
    """Generates pink noise using the Voss-McCartney algorithm.
    
    nrows: number of values to generate
    rcols: number of random sources to add
    
    returns: NumPy array
    """
    # Load the noise and define the mixture
    path = os.path.dirname(utilstf.__file__)
    babble,fs = librosa.load(os.path.join(path,'babble.wav'), sr=8000)
    Nb = len(babble)
    assert N<Nb, 'Cannot generate babble noise of length {}'.format(N)

    tinit = np.random.randint(0,Nb-N-1)

    return babble[tinit:tinit+N] - np.mean(babble[tinit:tinit+N])