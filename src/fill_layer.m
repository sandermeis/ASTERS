function layer = fill_layer(param, folderName)
arguments
    param
    folderName = "input"
end

fn = folderName + "/layers.xlsx";

opts = detectImportOptions(fn);
opts = setvartype(opts, 'input', 'string');  %or 'char' if you prefer
layer = table2struct(readtable(fn, opts, 'Sheet', param.layerSheet));

% Unsafe w.r.t. cybersecurity
run(folderName + "/" + param.surfaceFile + ".m")

% Loop through layers
for j = 1:numel(layer)
    layer(j).material = string(strsplit(layer(j).material,{' ',','}));
    layer(j).input = string(strsplit(layer(j).input,{' ',','}));

    % If input is numbers
    if all(ismember(char(layer(j).input), '123456789'))
        layer(j).input = num2cell(str2double(layer(j).input));
    % If input is uniform
    elseif layer(j).input=="u"||layer(j).input=="uniform"||layer(j).input=="Uniform"
        layer(j).input = 0;
    else
        error("Wrong layer input")
    end

    layer(j).add = param.add;
    layer(j).fill = param.fill;
    % Add Surfaces from the 'surfaces' cell array, based on the index entered in input
    if iscell(layer(j).input)
        for k = 1:numel(layer(j).input)
            layer(j).input{k} = surfaces{layer(j).input{k}};
        end
    end

end


%%% MAKE THEM WORK WITH STRING ARRAYS
layer = build_layerstack(layer, param);


% Check parameters
for i = 1:numel(layer)
    if iscell(layer(i).input)
        if param.res ~= size(layer(i).geometry.eps_struc, 1)
            error("Device resolution does not match input resolution")
        end
    end
end

% check
eps_lab = import_permittivities({layer.material}, param.wavelengthArray);
[layer.permittivities] = deal(eps_lab{:});
end
