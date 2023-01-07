function layer = fill_layer(param, folderName)
arguments
    param
    folderName = "input"
end

fn = folderName + "/layers.xlsx";

% Get layer from layers.xlsx and sheet parameter
opts = detectImportOptions(fn);
opts = setvartype(opts, 'input', 'string');
layer = table2struct(readtable(fn, opts, 'Sheet', param.layerSheet));

% Runs the code of the surface creation file
run(folderName + "/" + param.surfaceFile + ".m")

% Loop through layers and parse material and input
for j = 1:numel(layer)
    layer(j).material = string(strsplit(layer(j).material, {' ',','}));
    layer(j).input = string(strsplit(layer(j).input, {' ',','}));

    if ismissing(layer(j).input)
        error("Wrong layer input")
    % If input is numbers
    elseif all(ismember(char(layer(j).input), '123456789'))
        layer(j).input = num2cell(str2double(layer(j).input));
    % If input is uniform
    elseif layer(j).input == "u" || layer(j).input == "uniform" || layer(j).input == "Uniform"
        layer(j).input = 0;
    else
        error("Wrong layer input")
    end

    % Add Surfaces from the 'surfaces' cell array, based on the index entered in input
    if iscell(layer(j).input)
        if exist('surfaces', 'var')
            for k = 1:numel(layer(j).input)
                layer(j).input{k} = surfaces{layer(j).input{k}};
            end
        else
            error("Missing 'surfaces' cell array.")
        end
    end
end

% Creates final layers from surfaces input
layer = build_layerstack(layer, param);

% Check if surfaces match with parameters
for i = 1:numel(layer)
    if iscell(layer(i).input)
        if param.res ~= size(layer(i).geometry.eps_struc, 1)
            error("Device resolution does not match input resolution")
        end
    end
end

if param.calcFields
    total_thickness = sum([layer.L]);
    if total_thickness > 2500
        warning("Large total device thickness (" + string(total_thickness) + " nm) ,fields might not save correctly")
    end
end


% Get permittivities of the materials used and store in layer struct
eps_lab = import_permittivities({layer.material}, param.wavelengthArray);
[layer.permittivities] = deal(eps_lab{:});

end
