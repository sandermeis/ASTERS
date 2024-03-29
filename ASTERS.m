%%%%%%%%%% ASTERS %%%%%%%%%
% Simulation and process options can be found in the user guide.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear (Optional)
clear all
close all
clc

%% Add folders to path
addpath("src", "input", "src/objects", genpath("src/packages"), "src/refractive_indices", "src/data");

%% Load simulation parameters
param = load_parameters();

%% Run simulation
folderName = RCWA(param, simulationName="FullSolarCell");

%% Process results
batch = RCWA_process(folderName);