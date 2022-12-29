%%%%%%%%%% ASTERS %%%%%%%%%
% Simulation and process options can be found in the user guide.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear (Optional)
% clear all
% close all
% clc

%% Add folders to path
addpath("src/side_projects", "src", "src/data_daan", "src/data_maarten", "src/data_tom", "src/side_projects/generated_surfaces", "input", "src/objects", genpath("src/packages"), "src/refractive_indices");

%% Load simulation parameters
param = load_parameters();

%% Run simulation
folderName = RCWA(param);

%% Process results
RCWA_process(folderName)