# %%
import sys
sys.path.append("..")
import numpy as np
from matplotlib import pyplot as plt
from IPython.display import Audio
import scipy.signal as sg
import librosa
from src.utilities.utils import babble_noise
from mcsm_benchs.Benchmark import Benchmark
from pesq import pesq
from pystoi import stoi
import os
# from src.utilities.scale_functions import scale_fun_APF, scale_fun_Fvs, scale_fun_F

# from src.aps_metric.perf_metrics import musical_noise_measure_aps
from src.methods.method_hard_threshold import NewMethod as hard_thresholding
from src.methods.method_garrote_threshold import NewMethod as garrote_thresholding
from src.methods.method_delaunay_triangulation import delaunay_triangulation_denoising 

def dt_method(signal,*args,**kwargs):
    output = delaunay_triangulation_denoising(signal=signal,
                                            mcsim_path='src/utilities/simulations',
                                            *args,**kwargs)
    return output


# Create a dictionary of signals.
signals_dict = {}
N = 1*8192 # A bit more than a second of signal.
for i in range(6,7):
    x,fs = librosa.load('./signals/{}_male.mp3'.format(i), sr=8000)
    x = x[:N]
    signals_dict['{}_male'.format(i)]=x

for i in range(6,7):
    x,fs = librosa.load('./signals/{}_female.mp3'.format(i), sr=8000)
    x = x[:N]
    signals_dict['{}_female'.format(i)]=x

x,fs = librosa.load('./signals/cello.wav', sr=8000)
x = x[500:500+N]
signals_dict['cello']=x

x, fs = librosa.load(librosa.ex('trumpet'), duration=1.5, sr=8000)
x = x[0:N]
signals_dict['trumpet']=x


# %%
# Parameters
Nfft = 2*N
thr = np.arange(0.5,5.5,0.5)
lmax = np.arange(1.0,3.25,0.125)
SNRs = [-5, 0, 5, 10, 15, 20]
reps = 50

# %%
dict_methods = {
                'dt': dt_method,
                # 'ht': hard_thresholding().method,
                # 'st': garrote_thresholding().method,
                }

dict_parameters = {'dt': [{'LB':q} for q in lmax],}
# dict_parameters = {
                    # 'dt': [
                        #  {'scale_fun':scale_fun_APF},
                        #    {'scale_fun':scale_fun_Fvs},
                        #    {'scale_fun':scale_fun_F},
                        #    {'LB':q} for q in lmax,
                        #    ],
                    # 'ht': [{'coeff':q, 'Nfft':2**12} for q in thr],
                    # 'st': [{'coeff':q, 'Nfft':2**12} for q in thr],
                    # }

# %% Once new methods are added, load previous benchmark and add new methods.
new_methods_flag = False
# benchmark = Benchmark.load_benchmark('./results/benchmark_babble')
# benchmark.add_new_method(methods=dict_methods,parameters=dict_parameters)

# %% Here define performance functions.
# PESQ:
def pesq_metric(x,x_hat,**kwargs):
    return pesq(fs,x,x_hat,'nb')

# STOI:
def stoi_metric(x,x_hat,**kwargs):
    stoival = stoi(x, x_hat, fs, extended=False)
    return 100/(1+np.exp(-17.4906*stoival+9.6921))

from src.aps_metric.perf_metrics import musical_noise_measure_aps_8
def aps(x,xhat,**kwargs):
    return musical_noise_measure_aps_8(x,xhat,fs=fs,**kwargs)

print(stoi_metric(x,x), pesq_metric(x,x))

perf_funs = {'pesq':pesq_metric,
             'stoi':stoi_metric,
	        #  'aps':aps, # Only compute this at best pesq.
             }
if new_methods_flag:
    benchmark.objectiveFunction = perf_funs
    benchmark.verbosity = 5
    # benchmark.parallel_flag = 2

# %%
if not new_methods_flag:
    benchmark = Benchmark(task = 'denoising',
                            methods = dict_methods,
                            N = len(x), 
                            SNRin = SNRs, 
                            repetitions = 100, #reps,
                            signal_ids= signals_dict,
                            parameters=dict_parameters,
                            verbosity=5, 
                            parallelize=9,
                            obj_fun=perf_funs,
                            complex_noise=babble_noise,                     
                            )
# %%
benchmark.run() # Run the benchmark
benchmark.save_to_file(filename = './results/benchmark_babble_APF2')

# This formats the results on a DataFrame
results_parameters = benchmark.get_results_as_df()
results_parameters


