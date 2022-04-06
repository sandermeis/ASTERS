height = param.height;
numrandom = param.numrandom;
surfaceindex = param.surfaceindex;
featres = param.featres;
res = param.res;
size = param.size;

surfacetypes=["Cone","Pyramid","Sphere"];
surfacetype = surfacetypes(surfaceindex);

a = Surface(res, size);
f = Feature(featres, 2500, height, surfacetype);
a.addRandomFeatures(f, numrandom, "seed", 0, "PBC", true);
a.placeFeatures("PBC", true, "mode", "add");

b = Surface(256, 5000);
f2 = Feature(192, 2500, 20, "Sphere");
b.addRandomFeatures(f2, 12, "seed", 1, "PBC", true);
b.placeFeatures("PBC", true, "mode", "merge");

surfaces{1} = a;
surfaces{2} = b;