% TODO
%
% resolution matching feature/surface
% check perm outside range wavelengths
% Auto break on conv check
% Grating test
% image field inside layer
% check my saves on matlab file exchange
% maarten stuff, layers
% fix haze 4d figure

% TODO input
% Feature inputs range
% Fix useSurfaceSize
% Make work (give warning) with random strings after + that arent in params
% with just 1 paramset, ndgrid doesnt work
% repeat variable n times
% multiple wavelength ranges input
% vary both together (diagonal)
% "c" plus multiple arrays together

clear all
close all
%%
options.save = true;
options.parallel = true;
options.simulationName = 'MaartenTry1';

%%
[param, Sz] = RCWA(options);
%%
%RCWA_process(folderName)

