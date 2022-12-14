% To do list
%
% Error when harmonics and xy resolution are close, p0+pfft=0 as index
% wrong parameters when using +\param, param
% resolution matching feature/surface
% check perm outside range wavelengths
% Auto break on conv check
% Grating test
% make addfeature also accept list of obj
% Fix useSurfaceSize
% Make it so it copies whatever you named your 'createsurface', now just works for param(1)

%% Optional clearing
% clear all
% close all
% clc

%% Add file paths
addpath("src/side_projects","src","src/data_daan","src/data_maarten","src/data_tom","src/side_projects/generated_surfaces","input","src/objects",genpath("src/packages"),"src/refractive_indices");

%% Load simulation parameters
param = load_parameters();

%% Simulation options
options.simulationName = 'Simulation';
options.parallel = false;
options.onlinesave = false;
options.onlinePathName = '/run/user/1000/gvfs/smb-share:server=amsbackup-srv.science.ru.nl,share=amsbackup/Students/Sander/results/';

%% Run simulation
folderName = RCWA(param, options);

%% Process results
% Can input an array of foldernames
RCWA_process(folderName)