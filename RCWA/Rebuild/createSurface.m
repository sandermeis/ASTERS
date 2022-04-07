% height = param.height;
% numrandom = param.numrandom;
% surfaceindex = param.surfaceindex;
% featres = param.featres;
res = param.res;
size = param.size;

% surfacetypes=["Cone","Pyramid","Sphere"];
% surfacetype = surfacetypes(surfaceindex);

surf_algaas = Surface(res, size);
surf_oxide = Surface(res, size);
surf_oxide_uniform = Surface(res, size);

surf_algaas.addFeature(Feature("maarten/AlGaAs_surface1_10um.csv",10000),1,1);
surf_oxide.addFeature(Feature("maarten/Oxide_surface1_10um.csv",10000),1,1);
surf_oxide_uniform.addUniform(70);

% a = Surface(res, size);
% f = Feature(featres, 2500, height, surfacetype);
% a.addRandomFeatures(f, numrandom, "seed", 0, "PBC", true);
surf_algaas.placeFeatures("PBC", true, "mode", "add");
surf_oxide.placeFeatures("PBC", true, "mode", "add");
surf_oxide_uniform.placeFeatures("PBC", true, "mode", "add");

% b = Surface(256, 5000);
% f2 = Feature(192, 2500, 20, "Sphere");
% b.addRandomFeatures(f2, 12, "seed", 1, "PBC", true);
% b.placeFeatures("PBC", true, "mode", "merge");

surfaces{1} = surf_algaas;
surfaces{2} = surf_oxide;
surfaces{3} = surf_oxide;