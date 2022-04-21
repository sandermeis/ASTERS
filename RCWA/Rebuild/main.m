% TODO
%
% resolution matching feature/surface
% check perm outside range wavelengths
% Auto break on conv check
% Grating test
% image field inside layer
% check my saves on matlab file exchange
% maarten stuff, layers
% make addfeature also accept list of obj
% Fix useSurfaceSize

clear all
close all
%%
param = load_parameters();
%%
options.desktop = true;
options.onlinesave = true;
options.parallel = false;
options.simulationName = 'OxideUrbach25';

%%
folderName = RCWA(param, options);
%%
RCWA_process(folderName)

