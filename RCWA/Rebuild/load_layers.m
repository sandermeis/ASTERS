function [layerset, lay_ind] = load_layers(a,which_lay)

[l_C, ~, l_ic] = unique(which_lay);

lay_ind = num2cell(l_ic);

if ~isempty(l_C)
    fn = "layers.xlsx";
    %sheets = sheetnames(fn);
    for i = 1:numel(l_C)
        layerset(i).layer = table2struct(readtable(fn,'Sheet',l_C(i)));
    end

    for i = 1:numel(layerset)
        for j=1:numel(layerset(i).layer)
            if layerset(i).layer(j).input
                k = layerset(i).layer(j).input;
                if isnumeric(k) && k<=length(a) && k>=1
                    layerset(i).layer(j).input=a{k};
                else
                    error("Wrong layer input")
                end
            end
        end

        layerset(i).layer = build_layerstack(layerset(i).layer);
        eps_lab = import_permittivities({layerset(i).layer(:).material});
        [layerset(i).layer.permittivities] = deal(eps_lab{:});
    end
else
    error("No layers selected")
end

end

