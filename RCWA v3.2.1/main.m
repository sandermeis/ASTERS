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
% Make it so it copies whatever you named your 'createsurface', now just works for param(1)

% bb=readmatrix("refractive_indices/GaAs.csv");
% plot(2*pi./(bb(:,1)/1000), bb(:,1)/1000);
% hold on;
% plot(2*pi./(bb(:,1)/1000).*bb(:,2), bb(:,1)/1000)
% ylim([0.5,1])

%%

clear all
%close all
%%
addpath("data_daan","data_maarten","data_tom","side_projects","side_projects/generated_surfaces","input","objects",genpath("packages"),"refractive_indices");
%%
param = load_parameters();
%%
options.desktop = true;
options.onlinesave = false;
options.parallel = false;
options.simulationName = 'RetryBugFix';

%%
folderName = RCWA(param, options);
%%

RCWA_process(folderName)

% folderName = "MaartenFlatOxideTest2_24_05_22_11_38_13";
% RCWA_process(folderName)
% folderName = "MaartenTexturedOxideTest_24_05_22_12_10_13";
% RCWA_process(folderName)

% folderName = "SpacerThicknessTest_18_05_22_12_54_16";


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

