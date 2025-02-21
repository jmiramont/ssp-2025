""" Scale functions to compute the automatic scale in DT determination.
"""
import numpy as np
from src.utilities.spatstats_utils import compute_rank_envelope_test, generate_white_noise_zeros_pp
import pickle
from src.utilities.utils import apf_mc_test
import os

def scale_fun_Fvs(signal, mcsim_path=None):
    test_params = { 'alpha':0.05,
                        'fun':'Fest', 
                        'correction':'rs', 
                        'transform':'asin(sqrt(.))',
                        'rmin':0.65,
                        'rmax':1.05,                            
                    }
    scale_pp = compute_scale(signal, path=mcsim_path,**test_params)
    return scale_pp

def scale_fun_F(signal, mcsim_path=None):
    test_params = { 'alpha':0.05,
                        'fun':'Fest', 
                        'correction':'rs', 
                        # 'transform':'asin(sqrt(.))',
                        'rmin':0.65,
                        'rmax':1.05,                            
                    }
    scale_pp = compute_scale(signal, path=mcsim_path,**test_params)
    return scale_pp

def scale_fun_APF(signal, mcsim_path=None):
    output_dict = apf_mc_test(signal, hdim=1, return_all=True, nsim=199, path=mcsim_path)
    apf_signal = output_dict['apf_obs']
    apf_noise = output_dict['apf_noise']
    # apf_mean = output_dict['apf_mean']
    ms = output_dict['ms']
    # nsims = output_dict['nsims']
    # upper_env = np.max(apf_noise,axis=0)
    lower_env = np.min(apf_noise,axis=0)    
    # r_max_dif = np.max(lower_env-apf_signal)
    ind_max_dif = np.argmax(lower_env-apf_signal)
    ms2 = np.sqrt(ms)
    scale_pp = ms2[ind_max_dif]
    return scale_pp

def compute_scale(signal, path='', **test_params):
    """_summary_

    Args:
        signal (_type_): _description_
        Nfft (_type_): _description_
        cs (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    
    # If the simulated ppp are present, use them, otherwise, generate the simulation and
    # save it so that it can be used later (this saves time, avoiding re-simulation).
    nsim = 2499
    N = len(signal)
    try:
        with open(os.path.join(path,'ppp_simulations_N_{}_Nsim_{}.mcsim'.format(N,nsim)), 'rb') as handle:
            ppp_simulation = pickle.load(handle)
    except:
        print('MC Simulation running...')
        list_ppp = generate_white_noise_zeros_pp(N, nsim=nsim)
        ppp_simulation = {'list_ppp':list_ppp,
                    'N' : N,
                    'nsim' : nsim
                    }
        with open(os.path.join(path,'ppp_simulations_N_{}_Nsim_{}.mcsim'.format(N,nsim)), 'wb') as handle:
            pickle.dump(ppp_simulation, handle, protocol=pickle.HIGHEST_PROTOCOL)   


    output_dic = compute_rank_envelope_test(signal,
                                         return_dic=True, 
                                         ppp_sim=ppp_simulation['list_ppp'],
                                         **test_params)

    reject_H0 = output_dic['rejectH0']
    if not reject_H0:
        print('Radius computed but signal was not detected.')
    # if not reject_H0:
    #     # print('No detection.')
    #     radius_of_rejection = 0.8
    # else:
    radius_of_rejection = output_dic['r_max_dif']
    
    return radius_of_rejection