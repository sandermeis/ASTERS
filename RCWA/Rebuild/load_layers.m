function [layer_structure, lay_ind, pset] = load_layers(p, which_lay)

num_p = numel(p);

surfaces = cell(1, num_p);
% unique combination of params
for i = 1:num_p

    surfacetypes=["Cone","Pyramid","Sphere"];

    height = p(i).p1;
    numrandom = p(i).p2;
    surfaceindex = p(i).p3;
    res = p(i).p4;

    surfacetype = surfacetypes(surfaceindex);

    a = Surface(512, 10000);
    f = Feature(res, 5000, height, surfacetype);

    a.addRandomFeatures(f, numrandom, "seed", 0, "PBC", true); % add option for random seed

    a.placeFeatures("PBC", true, "mode", "merge");

    surfaces{i} = a;

end

pset = num2cell(1:num_p);
% combine all into surf
% pset becomes index of generated surf (input)

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Which layer stacks "which sheets in excel"
[l_C, ~, l_ic] = unique(which_lay);

lay_ind = num2cell(l_ic);

if ~isempty(l_C)
    fn = "layers.xlsx";
    %sheets = sheetnames(fn);
    for i = 1:numel(l_C)
        for j = 1:num_p
            layer_structure(i).layer_set(j).layer = table2struct(readtable(fn,'Sheet',l_C(i)));
        end
    end

    for i = 1:numel(layer_structure)
        % Which layers "surfaces"
        % now: make new layer for each new parameter
        %         for j=1:numel(layer_structure(i).layer)
        %             if layer_structure(i).layer(j).input
        %                 k = layer_structure(i).layer(j).input;
        %                 if isnumeric(k) && k<=length(a) && k>=1
        %                     layer_structure(i).layer(j).input=a{k};
        %                 else
        %                     error("Wrong layer input")
        %                 end
        %             end
        %         end

        for j=1:num_p
            for k=1:numel(layer_structure(i).layer_set(j).layer)
                if layer_structure(i).layer_set(j).layer(k).input
                    layer_structure(i).layer_set(j).layer(k).input = surfaces{j};
                elseif ~layer_structure(i).layer_set(j).layer(k).input

                else
                    error("Wrong layer input")
                end
            end


        layer_structure(i).layer_set(j).layer = build_layerstack(layer_structure(i).layer_set(j).layer);
        eps_lab = import_permittivities({layer_structure(i).layer_set(j).layer(:).material});
        [layer_structure(i).layer_set(j).layer.permittivities] = deal(eps_lab{:});
        end
    end
else
    error("No layers selected")
end

end

