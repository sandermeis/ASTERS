% TODO
% 
% height distribution
% PSD, ACF
% Surface remove X and Y for RAM saving
% uiconfirm
% parallel/save options
% 
%
% resolution matching feature/surface
% Auto break on conv check
% Grating test
% image field inside layer
% check my saves on matlab file exchange
% maarten stuff, layers
% fix haze 4d figure

% TODO input
% Make work (give warning) with random strings after + that arent in params
% with just 1 paramset, ndgrid doesnt work
% repeat variable n times
% multiple wavelength ranges input
% vary both together (diagonal)
% "c" plus multiple arrays together

clear all
close all
%%
% save yes/no, display graphs, yes/no
options.save = true;
options.parallel = true;
%%
folderName = RCWA(options);
%%
%RCWA_process(folderName)

