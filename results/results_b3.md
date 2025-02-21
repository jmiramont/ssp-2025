# Benchmark Report [[Results .csv]](https://jmiramont.github.io/ssp-2025/results/b3/results.csv) 


## Configuration

Length of signals: 8192

Repetitions: 50

SNRin values: 
10, 


### Methods  

* dt 

* ht 

### Signals  

* 6_female 

## Mean results tables: 

The results shown here are the average and 95\% CI of                               the performance metric with Bonferroni correction.                               Best performances are **bolded**. 
### Signal: 6_female [[View Plot]](https://jmiramont.github.io/ssp-2025/results/b3/plot_6_female.html)    [[Get .csv]](https://jmiramont.github.io/ssp-2025/results/b3/results_6_female.csv)
|    | Method + Param                                              | SNRin=10dB (average)   | SNRin=10dB (CI)    |
|---:|:------------------------------------------------------------|:-----------------------|:-------------------|
|  0 | dt{'scale_fun': <function scale_fun_APF at 0x7d29c335fe20>} | 40.30                  | ['35.08', '44.99'] |
|  1 | dt{'scale_fun': <function scale_fun_Fvs at 0x7d29c43dea70>} | 61.76                  | ['55.04', '67.39'] |
|  2 | dt{'scale_fun': <function scale_fun_F at 0x7d29c335fd90>}   | **67.15**              | ['63.86', '70.12'] |
|  3 | ht                                                          | 64.74                  | ['59.39', '69.62'] |
