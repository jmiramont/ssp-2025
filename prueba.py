# %% [markdown]
# # Displaying results for the signal detection benchmark
# 
# This notebooks loads the results of the benchmark and generates the figures shown in the paper.

# %%
from mcsm_benchs.Benchmark import Benchmark
from mcsm_benchs.ResultsInterpreter import ResultsInterpreter
from mcsm_benchs.SignalBank import SignalBank

import scipy.stats as spst
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from numpy import pi as pi
import scipy.signal as sg
import matplotlib.pyplot as plt
from src.utilities.utilstf import get_spectrogram
import os


# %%
# Load the benchmark results
filename = os.path.join('results','last_benchmark_detection_1024_new')
benchmark1 = Benchmark.load_benchmark(filename)

filename = os.path.join('results','last_benchmark_detection_APF')
benchmark2 = Benchmark.load_benchmark(filename)


# %%
benchmark2.objectiveFunction

# %%

benchmark = benchmark1+benchmark2
# benchmark.name = 'Single Chirp Detection - N=1024'
# benchmark.description = "Detecting a chirp in white Gaussian noise."
# benchmark.save_to_file(filename)

print(benchmark.name, benchmark.description)


# %%
# Load interpreter.
interpreter = ResultsInterpreter(benchmark)
df = interpreter.rearrange_data_frame()
print(np.unique(df['Method']))

# %%
# Check DataFrame of results.
# df = interpreter.get_benchmark_as_data_frame()
df = interpreter.rearrange_data_frame()

df = df[df['Parameter']!="{'statistic': 'Frs', 'pnorm': 2, 'rmax': 0.5, 'MC_reps': 2499}"]
df = df[df['Parameter']!="{'statistic': 'Frs_vs', 'pnorm': 2, 'rmax': 0.5, 'MC_reps': 2499}"]
df = df[df['Parameter']!="{'statistic': 'Frs_vs', 'pnorm': 2, 'rmax': 2.0, 'MC_reps': 2499}"]
df = df[df['Parameter']!="{'statistic': 'Frs', 'pnorm': 2, 'rmax': 2.0, 'MC_reps': 2499}"]
df = df[df['Parameter']!="{'statistic': 'Frs', 'pnorm': 2, 'rmax': 1.0, 'MC_reps': 199}"]
df = df[df['Parameter']!="{'statistic': 'Frs_vs', 'pnorm': 2, 'rmax': 1.0, 'MC_reps': 199}"]


df = df[df['Parameter']!="{'statistic': 'Frs', 'MC_reps': 199}"]
df = df[df['Parameter']!="{'statistic': 'Frs_vs', 'MC_reps': 199}"]

np.unique(df['Parameter'])

# %%
# Use this function only for the CP CI shown in the interactive figures using Plotly:
def clopper_pearson(x, alpha=0.05, bonferroni=1):
    """
    Clopper-Pearson confidence interval for Bernoulli parameter
    alpha: confidence level
    k: number of successes
    n: number of observations
    """
    alpha = alpha/bonferroni
    n = len(x) # k: number of successes
    k = sum(x) # n: number of observations
    lb = np.mean(x) - spst.beta.ppf(alpha/2, k, n-k+1)
    ub = spst.beta.ppf(1 - alpha/2, k+1, n-k)-np.mean(x)
    return lb, ub

cp_ci = lambda x: clopper_pearson(x, alpha=0.05, bonferroni=3)

# %%
# Report shown in the repo 
interpreter.save_report(path='../results', 
                        link='https://jmiramont.github.io/benchmarks-detection-denoising/results/detection')

# Interactive figures shown in the repo
interpreter.get_html_figures(df=interpreter.get_benchmark_as_data_frame(), varfun=cp_ci, path='../results/detection', bars=True, ylabel='Detection Power')

# .csv files for sharing results
interpreter.get_csv_files(path='../results/detection')

# %%
from plotly.offline import iplot
fig = interpreter.get_summary_plotlys(varfun=cp_ci)

for f in fig:
    f.update_layout(yaxis_title='Power', xaxis_title='SNRin (dB)')
    iplot(f)

# %%
# Use this function for the CP CI shown in the interactive figures using matplotlib:
def clopper_pearson(x, alpha=0.05, bonferroni=1):
    """
    Clopper-Pearson confidence interval for Bernoulli parameter
    alpha: confidence level
    k: number of successes
    n: number of observations
    """
    alpha = alpha/bonferroni
    n = len(x) # k: number of successes
    k = sum(x) # n: number of observations
    lb = spst.beta.ppf(alpha/2, k, n-k+1) 
    ub = spst.beta.ppf(1 - alpha/2, k+1, n-k)
    return lb, ub

cp_ci = lambda x: clopper_pearson(x, alpha=0.05, bonferroni=8)

# %%
# Fix new legends.
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

new_legends = [
            r'Global MAD-$F$',
            r'Global MAD-$\widetilde{F}$',
            r'Rank-$\widetilde{F}$-$r\in[0.65;1.05]$',
            r'Rank--$F$-$r\in[0.65;1.05]$',
            r'Rank-$\widetilde{F}$',
            r'Rank--$F$',
            r'MC-$F$-$r_{\max}=1.0$',
            r'MC-$\widetilde{F}$-$r_{\max}=1.0$',
            ]

# Set figure size (in inches):
fig_size_w = 8.2
fig_size_h = 2

errbar_params = {'errwidth':0.5,
                'capsize':0.3,
                'gap':0.1
                }

figs = interpreter.get_summary_plots(   df_rearr=df,
                                        size = (fig_size_w,fig_size_h),
                                        # filter_crit='any',
                                        # filter_str=['global',"'pnorm': 2, 'rmax': 1.0,",],
                                        errbar_fun=cp_ci, 
                                        savetofile=False, 
                                        plot_type='bars',
                                        errbar_params=errbar_params
                                        )

for i, fig in enumerate(figs):
    # Get signal the signal for each figure and compute spectrogram
    ax = fig.axes[0]
    s = benchmark.signal_dic['LinearChirp']
    S, stft  = get_spectrogram(s)
    ax = fig.axes[0]

    # Set inset axis with the spectrogram of the signal
    axins = inset_axes(ax, width=0.85, height=0.85, loc=2)
    axins.imshow(S, origin='lower')
    # axins.axis('off')
    fig.canvas.draw()
    axins.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)

    #Uncomment to remove legends from figure (they are print in another figure)
    ax.get_legend().remove()
    
    # h, l  = ax.get_legend_handles_labels()
    # ax.legend(h, new_legends, fontsize='small', frameon=False, loc=(0, 1.1), ncol=4)
    ax.set_title('N={}'.format(interpreter.N),fontsize=8)
    ax.set_ylabel(ax.get_ylabel(),fontsize=8)
    ax.set_xlabel(ax.get_xlabel(),fontsize=8)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=7)
    ax.set_ylim((0,1.05))
    fig.canvas.draw()
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
    
    # Save figure
    filename = os.path.join('figures','power_'+str(i)+'_N_{}.pdf'.format(interpreter.N))
    fig.savefig(filename, dpi=900, transparent=False, bbox_inches='tight')

# Get legends in a different figure
legendFig = plt.figure()
legendFig.set_size_inches((fig_size_w,fig_size_h))

h,_ = figs[0].axes[0].get_legend_handles_labels()
legendFig.legend(h, new_legends, fontsize='small', frameon=False, loc='center', ncol=4)
legendFig.canvas.draw()

# Save figure with legends
filename = os.path.join('figures','legend_power.pdf')
legendFig.savefig(filename, dpi=900, transparent=False, bbox_inches='tight')

plt.show()


