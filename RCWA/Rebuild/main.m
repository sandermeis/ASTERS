% TODO
% 
% height distribution
% PSD, ACF
% Surface remove X and Y for RAM saving
% uiconfirm
% save/parallel
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
param = load_parameters();

%%
p = struct('p1', {param.p1}, 'p2', {param.p2},'p3', {param.p3},'p4', {param.p4});
[layer_structure, lay_ind, pset] = load_layers(p,[param.lay]);

[param.pset] = pset{:};
[param.lay] = lay_ind{:};
%%
% save yes/no, display graphs, yes/no
options.save = true;
options.parallel = true;
%%
Sz = RCWA(layer_structure, param, options);

