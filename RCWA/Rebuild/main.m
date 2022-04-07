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
% copy createsurface file into results
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
options.save = false;
options.parallel = false;
options.simulationName = 'MaartenOxidesTextured';

%%
[param, Sz] = RCWA(options);
%%
RCWA_process(Sz)

