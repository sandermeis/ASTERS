% TODO
%
% Error when harmonics and xy resolution are close, p0+pfft=0 as index
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
addpath("data_daan","data_maarten","input","objects",genpath("packages"),"refractive_indices");
%%
param = load_parameters();
%%
options.desktop = true;
options.onlinesave = false;
options.parallel = true;
options.simulationName = 'WedgeAlGaAs2kConvergenceTest';

%%
folderName = RCWA(param, options);
%%
folderName = "SpacerThicknessTest_18_05_22_12_54_16";
RCWA_process(folderName)

