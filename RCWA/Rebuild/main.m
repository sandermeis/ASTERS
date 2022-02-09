clearvars
close all

%%% SHAPE
%%%
%%% 0 = Uniform
%%% 1 = GratingX
%%% 2 = GratingY
%%% 3 = GratingXY
%%% 4 = Triangle
%%% 5 = Circle
%%% 6 = Rough
%%% 7 = Real surface, only 512x512, 10x10 or 30umx30um

% layer               = cell2struct(...
%                         [{'MgF2', 'ZnS', 'AlInP', 'GaAs', 'InGaP', 'Al03GaAs', 'Ag'};...
%                         {94, 44, 20, 300, 100, 100, 100};...
%                         {0, 0, 0, 0, 0, 0, 0};...
%                         {0, 0, 0, 0, 0, 0, 0}], {'material','L','shape','roughdim'}, 1);
% 
layer               = cell2struct(...
                        [{'InGaP', 'Al03GaAs', 'InGaP', 'GaP', 'Ag', 'Ag'};...
                        {125, 50, 50, 360, 213, 100};...
                        {0, 0, 0, 0, 6, 0};...
                        {0, 0, 0, 0, 20, 0}], {'material','L','shape','roughdim'}, 1);

% layer               = cell2struct(...
%                         [{'GaAs','GaAs'};...
%                         {100,500};...
%                         {1, 0};...
%                         {0, 0}], {'material','L','shape','roughdim'}, 1);

wavelengthArray     = 650:5:700;

layer               = import_permittivities(layer,wavelengthArray,false);

input_wave.theta	= 0;%8/360*2*pi;
input_wave.phi      = 0;
input_wave.pTE      = 0.5;
input_wave.pTM      = 0.5;

device.ref_medium	= "Air";
device.trn_medium	= "Air";
device.eps_ref      = 1;
device.mu_ref       = 1;
device.eps_trn      = 1;
device.mu_trn       = 1;
device.size_x       = 30000;
device.size_y       = 30000;
device.res_x        = 512;
device.res_y        = 512;

device.k_0          = 2*pi./wavelengthArray;
device.truncateHarm = true;
device.truncFactor  = 1;
device.calcAllRough = false;
device.checkConverg = false;
device.Hmax         = 5;

% TODO
% 1. Negative absorption rough layer
% 2. Improve visualisation
% 3. Checking convergence
% 4. Averaging
% 5. Improve truncation
% 6. Creating shapes
% 7. Anechoic chamber shape
% 8. Grating test
% image field inside layer
% check my saves on matlab file exchange

[device, layer, input_wave, Sz] = RCWA(layer,device,input_wave,wavelengthArray);
%%
RCWA_plot(layer, device, wavelengthArray, Sz);