% Create surfaces here, optionally based on the used defined parameters in
% input.txt, found in the struct "param"
% -------------------------------------------------------------------------

% surf_algaas = Surface(param.res, param.size);
% surf_oxide = Surface(param.res, param.size);
% surf_oxide_uniform = Surface(param.res, param.size);
% 
% surf_algaas.addFeature(Feature("data_maarten/AlGaAs_surface1_10um.csv",10000),1,1);
% surf_oxide.addFeature(Feature("data_maarten/Oxide_surface1_10um.csv",10000),1,1);
% surf_oxide_uniform.addUniform(70);
% 
% surf_algaas.placeFeatures("PBC", true, "mode", "add");
% surf_oxide.placeFeatures("PBC", true, "mode", "add");
% surf_oxide_uniform.placeFeatures("PBC", true, "mode", "add");

% Add surfaces here, index can be used in the "input" field in layers.xlsx
% -------------------------------------------------------------------------

% surfaces{1} = surf_algaas;
% surfaces{2} = surf_oxide;
% surfaces{3} = surf_oxide_uniform;

% Modify "layer" struct here
% -------------------------------------------------------------------------
layer(3).material = "GaAs_"+string(param.urbach)+"meV";
layer(4).L = param.algaasthickness;
layer(5).L = param.oxidethickness;


