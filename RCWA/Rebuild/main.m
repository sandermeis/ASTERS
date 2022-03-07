clear all
close all
%%
% TODO
%
% overarching parameter space prepare, with saving data
% check 1D features
% Averaging
% resolution matching
% More shapes
% resize with interpolation
% Grating test
% report RMS roughness
% image field inside layer
% check my saves on matlab file exchange
% save data, to average later manually

%%% Layer shape
%%%
%%% Features:
%%% "GratingX"
%%% "GratingY"
%%% "GratingXY"
%%% "Triangle"
%%% "Circle" 2D Circle, 3D Cylinder
%%% "Sphere" 3D Half sphere
%%% "Pyramid" 3D Pyramid
%%% "RidgeX" 3D Ridge X
%%% "RidgeY" 3D Ridge Y
%%% "WedgeX" 3D Wedge X
%%% "WedgeY" 3D Wedge Y
%%% "Cone" 3D Cone
%%% "Custom" Custom .csv input

%%% Surface:
%%% Rough1
%%% Rough2

%     650
%         Z = RD{3,1}; %650 %max 487
%     700
%         Z = RD{3,4}; %700 %max 1800, intermax 700
%     730
%         Z = RD{3,7}; %730 %max 213

% FSC
% layer               = cell2struct(...
%                         [{'MgF2', 'ZnS', 'AlInP', 'GaAs', 'InGaP', 'Al03GaAs', 'Ag'};...
%                         {94, 44, 20, 300, 100, 100, 100};...
%                         {0, 0, 0, 0, 0, 0, 0};...
%                         {0, 0, 0, 0, 0, 0, 0}], {'material','L','shape','roughdim'}, 1);

% Lara
% layer               = cell2struct(...
%                         [{'MgF2', 'ZnS', 'AlInP', 'GaAs', 'InGaP', 'MgF2', 'Ag'};...
%                         {94, 44, 25, 2100, 100, 100, 3000};...
%                         {0, 0, 0, 0, 0, 0, 0};...
%                         {0, 0, 0, 0, 0, 0, 0}], {'material','L','shape','roughdim'}, 1);

% Daan
% layer               = cell2struct(...
%                         [{'InGaP', 'Al03GaAs', 'InGaP', 'GaP', 'Ag', 'Ag'};...
%                         {125, 50, 50, 100, 213, 100};...
%                         {0, 0, 0, 0, 7, 0};...
%                         {0, 0, 0, 0, 50, 0}], {'material','L','shape','roughdim'}, 1);

% randsurf = Surface(10000,512);
% randsurf.addRoughsurf('mode','Artificial','hurst',1)
% randsurf.listFeatures()
% randsurf.placeFeatures()
% randsurf.plot()

% Maarten
% surf_algaas=Surface(10000,512);
% surf_oxide=Surface(10000,512);
% 
% surf_algaas.addFeature(Feature("maarten/AlGaAs_surface1_10um.csv",10000),1,1);
% surf_oxide.addFeature(Feature("maarten/Oxide_surface1_10um.csv",10000),1,1);
% 
% surf_algaas.placeFeatures()
% surf_oxide.placeFeatures()
% 
% surf_algaas.listFeatures()
% surf_oxide.listFeatures()
% 
% % surf_algaas.plot();
% % surf_oxide.plot();


for n=1:40

a=Surface(10000,512);
f_titan = Feature(round(512/2),5000,500,"Cone");
f_giant = Feature(round(512/2.5),5000,450,"Cone");
f_large = Feature(round(512/3.5),2500,400,"Cone");
f_medium = Feature(round(512/5),1250,300,"Cone");
f_small = Feature(round(512/10),625,100,"Cone");

a.addRandomFeatures(f_titan,2,"PBC",true)
a.addRandomFeatures(f_giant,5,"PBC",true)
a.addRandomFeatures(f_large,15,"PBC",true)
a.addRandomFeatures(f_medium,80,"PBC",true)
a.addRandomFeatures(f_small,100,"PBC",true)

a.placeFeatures("PBC",true, "mode", "merge");
a.listFeatures()
a.plot()
% b4 = Feature("Feature1_2500nm.csv",2500);
% 
% b5 = Feature("Feature1_2500nm.csv",2500);
% b5 = rotate(b5);
%a.addFeature(b,1,1);
% a.addFeature(b5,1,1)
% a.addFeature(b2,1,1)
% a.addFeature(b3,1,1)
% 

% a.placeRandomFeatures(1,1000,"PBC",true, "mode", "merge")
%a.hscale(200)
%a.addRandomFeatures(b4,3,"PBC",true)
%a.addRandomFeatures(b5,2,"PBC",true)
% a.placeFeatures("PBC",true, "mode", "merge");
%a.placeRandomFeatures(2,3,"PBC",true, "mode", "merge")
%a.plot();
%%

% Maarten

% layer               = cell2struct(...
%                         [{'GaAs','Al03GaAs','InGaP','GaAs'};...
%                         {100,200,100,5000};...
%                         {a,0, 0, 0};...
%                         {10, 0, 0, 0}], {'material','L','input','roughdim'}, 1);

% Test
% layer               = cell2struct(...
%                         [{'GaAs','GaAs','Ag'};...
%                         {500,200,1000};...
%                         {a,0, 0};...
%                         {20,0, 0};...
%                         {7,0, 0}], {'material','L','input','roughdim','shape'}, 1);

% Daan paper
layer               = cell2struct(...
                        [{'MgF2','ZnS','AlInP','GaAs','InGaP','Al03GaAs','InGaP','GaP','Ag'};...
                        {94,44,25, 300, 110, 130, 50, 500, 3000};...
                        {0, 0,0,0, 0, 0, 0, 0, a};...
                        {0, 0,0,0, 0, 0, 0, 0, 40}], {'material','L','input','roughdim'}, 1);


wavelengthArray     = 850:10:900;

layer               = import_permittivities(layer,wavelengthArray, 'plot_permfig', false);

% input wave
input_wave.theta	= 8/360*2*pi;
input_wave.phi      = 0;
input_wave.pTE      = 0.5;
input_wave.pTM      = 0.5;

% left and right media
device.ref_medium	= "Air";
device.trn_medium	= "Air";
device.eps_ref      = 1;
device.mu_ref       = 1;
device.eps_trn      = 1;
device.mu_trn       = 1;

% Resizes, requires interpolation
device.useSurfaceSize  = true;
device.size_x       = 10000;
device.size_y       = 10000;
device.res_x        = 512;
device.res_y        = 512;

device.k_0          = 2*pi./wavelengthArray;

% Truncate harmonics
device.truncateHarm = true;
device.truncFactor  = 1;
device.truncfig     = false;

% Calculate field in all rough layers
device.calcAllRough = false;

% Convergence check
device.checkConverg = false;
device.Hmax         = 9;

% Rough
device.plotSurf     = false;
device.reverse      = true;
device.recalcRoughL = true;

% Merge layers that are approximately equal
device.optimRough   = false;
device.tolerance    = 5e-3;


%%
[device, layer, input_wave, Sz] = RCWA(layer,device,input_wave,wavelengthArray);
filename = "overnight/test"+n;
save(filename)
clearvars -except n
close all
end
%%
%h = RCWA_plot(layer, device, wavelengthArray, Sz);

%%
%plot(wavelengthArray,squeeze(sum(Sz(:,8,:),1)))