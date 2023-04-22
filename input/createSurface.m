%% Create surfaces here
% Optionally based on the used defined parameters in input.txt, found in
% the struct "param"
% -------------------------------------------------------------------------


%% Add surfaces here
% Index can be used in the "input" field in layers.xlsx
% Example: Add two surfaces ("surf_oxide" and "surf_algaas") to the list,
% which can be accessed in layers.xlsx by entering 1 or 2 for surf_oxide
% and surf_algaas respectively.
% surfaces{1} = surf_oxide;
% surfaces{2} = surf_algaas;
% -------------------------------------------------------------------------


%% Modify "layer" struct here
% Example: Modify thickness of second layer based on the parameter "oxidethickness".
% layer(2).L = param.oxidethickness;
% -------------------------------------------------------------------------


%% Example 2
% surf_r = Surface(param.res, param.size);
% surf_r.addRoughsurf('height', 100);
% surf_r.placeFeatures();
% surf_r.plot
% surf_r.report
% surfaces{1} = surf_r;

%% Example 3
% surf_sph = Surface(param.res, param.size);
% f1 = Feature(param.res/2, param.size/2, 100, "Sphere");
% f2 = Feature(param.res/4, param.size/4, 50, "Sphere");
% f3 = Feature(param.res/8, param.size/8, 25, "Sphere");
% surf_sph.addRandomFeatures(f1, 5);
% surf_sph.addRandomFeatures(f2, 25);
% surf_sph.addRandomFeatures(f3, 125);
% surf_sph.placeFeatures();
% surf_sph.plot
% surf_sph.report
% surfaces{1} = surf_sph;

%% Example 4
% surf_sph = Surface(param.res, param.size);
% f1 = Feature(param.res/2, param.size/2, 100, "Sphere");
% f2 = Feature(param.res/4, param.size/4, 50, "Sphere");
% f3 = Feature(param.res/8, param.size/8, 25, "Sphere");
% surf_sph.addRandomFeatures(f1, 5, "seed", param.random_seed);
% surf_sph.addRandomFeatures(f2, 25, "seed", param.random_seed);
% surf_sph.addRandomFeatures(f3, 125, "seed", param.random_seed);
% surf_sph.placeFeatures();
% surf_sph.plot
% surf_sph.report
% surfaces{1} = surf_sph;

%% Example 5
% surf_r = Surface(param.res, param.size);
% surf_r.addRoughsurf('height', 100);
% surf_r.placeFeatures();
% surf_uniform = Surface(param.res, param.size);
% surf_uniform.addUniform(50);
% surf_uniform.placeFeatures();
% surfaces{1} = surf_r;
% surfaces{2} = surf_uniform;