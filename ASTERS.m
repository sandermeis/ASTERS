% To do list

% show fields input
% fix field resolution,
% load_parameters: warning only 1 wl is possible

% res must be at least 2x Hmax, MAKE A CHECK FOR THIS, option to duplicate array
% Error when harmonics and xy resolution are close, p0+pfft=0 as index
%% Optional clearing
clear all
close all
clc

%% Add file paths
addpath("src/side_projects","src","src/data_daan","src/data_maarten","src/data_tom","src/side_projects/generated_surfaces","input","src/objects",genpath("src/packages"),"src/refractive_indices");

%% Load simulation parameters
param = load_parameters();

%% Run simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation options can be found in the guide
% simulationName = 'Simulation';
% parallel = true;
% onlinesave = true;
% onlinePathName = 'YourPathName';
%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderName = RCWA(param);

%% Process results
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process options
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichsims = [];
% which_layer_jsc = "GaAs";
% which_plot = "All";
% disp_permittivity = [];
% disp_truncation = [];
% plot_layers = [];
% Can input an array of foldernames
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Average over sims

RCWA_process(folderName, whichsims = [], which_layer_jsc = "Hirst/GaAs")