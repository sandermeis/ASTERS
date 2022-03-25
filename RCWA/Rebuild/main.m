% TODO
% 

% save data, to average later manually
% Jsc in workspace
% height distribution
% PSD, ACF
% Surface remove X and Y for RAM saving
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

clear all
close all
%%
param = load_parameters();

%%
a1 = Surface(512,10000);
a2 = Surface(512,10000);
a3 = Surface(512,10000);

f_cone = Feature(512,5000,500,"Cone");
f_pyramid = Feature(512,5000,500,"Pyramid");
f_sphere = Feature(512,5000,500,"Sphere");

a1.addRandomFeatures(f_cone,10,"PBC",true) % add option for random seed
a2.addRandomFeatures(f_pyramid,10,"PBC",true)
a3.addRandomFeatures(f_sphere,10,"PBC",true)

a1.placeFeatures("PBC",true, "mode", "merge");
a2.placeFeatures("PBC",true, "mode", "merge");
a3.placeFeatures("PBC",true, "mode", "merge");

surfaces={a1,a2,a3};
%%
[layerset, lay_ind] = load_layers(surfaces,[param.lay]);
[param.lay] = lay_ind{:};
%%
% save yes/no, display graphs, yes/no
options.save = true;
options.parallel = true;
%%
Sz = RCWA(layerset, param, options);

