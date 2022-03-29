function [layer] = fill_layer(p, n, which_lay)

surfacetypes=["Cone","Pyramid","Sphere"];

height = p(n).p1;
numrandom = p(n).p2;
surfaceindex = p(n).p3;
res = p(n).p4;

surfacetype = surfacetypes(surfaceindex);

a = Surface(512, 10000);
f = Feature(res, 5000, height, surfacetype);

a.addRandomFeatures(f, numrandom, "seed", 0, "PBC", true); % add option for random seed

a.placeFeatures("PBC", true, "mode", "merge");

surfaces{1} = a;

% combine all into surf
% pset becomes index of generated surf (input)

%%%%%%%%%%%%%%%%%%%%%%%

% Which layer stacks "which sheets in excel"

fn = "layers.xlsx";

layer = table2struct(readtable(fn,'Sheet',which_lay));

for j=1:numel(layer)
    if layer(j).input
        k = layer(j).input;
        if isnumeric(k) && k<=length(a) && k>=1
            layer(j).input = surfaces{k};
        else
            error("Wrong layer input")
        end
    end
end

layer = build_layerstack(layer);
eps_lab = import_permittivities({layer.material});
[layer.permittivities] = deal(eps_lab{:});
end
