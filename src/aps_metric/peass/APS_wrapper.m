function [APS] = APS_wrapper(xi,xn,xo,fs)

options.destDir = './';
options.segmentationFactor = 1; % increase this integer if you experienced "out of memory" problems

% Save inputs as files for the PEASS function:
dir = './';

rng('shuffle')
cod = string(randi(15000));

xi_file = char('tmp_xi_'+cod+'.wav');
audiowrite([dir xi_file],xi/max(abs(xi)),fs);
path_in{1} = [dir xi_file];

if ~ isempty(xn)
    xn_file = char('tmp_xn_'+cod+'.wav');
    audiowrite([dir xn_file],xn/max(abs(xn)),fs);
    path_in{2} = [dir xn_file];
end

xo_file = char('tmp_xo_'+cod+'.wav');
audiowrite([dir xo_file],xo/max(abs(xo)),fs);
path_out = [dir xo_file];

res = PEASS_ObjectiveMeasure(path_in,path_out,options);
APS = res.APS;

system(char('rm tmp_*_'+cod+'*'));


