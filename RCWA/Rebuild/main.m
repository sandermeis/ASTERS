% TODO
%
% Error when harmonics and xy resolution are close, p0+pfft=0 as index
% wrong parameters when using +\param, param
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

% foldernames = ["DuplicateSurfaceTestSteppedGratingGaAsConvergence_18_05_22_14_45_33", "DuplicateSurfaceTestSteppedGratingGaAsTry2_18_05_22_14_01_44","DuplicateSurfaceTestSteppedGratingGaAsConvergence1113_18_05_22_15_45_18"];
% 
% RCWA_process(foldernames)
% folderName = "Presentation1_17_05_22_18_20_24";

% folderName = "DuplicateSurfaceTestDiffangle_18_05_22_11_50_02";
% RCWA_process(folderName)

% folderName = "Presentation1_17_05_22_18_20_24";
% RCWA_process(folderName)
% 
% folderName = "PresentationOxideSpacer2_18_05_22_09_19_59";
% RCWA_process(folderName)
% 
% folderName = "PresentationOxideSpacer_17_05_22_22_40_17";
% RCWA_process(folderName)

