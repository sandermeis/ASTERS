function RCWA_process(folderName, plot_options)
arguments
    folderName (1,:) string
    plot_options.whichsims (1,:) {mustBeInteger} = []
    plot_options.which_layer_jsc (1,1) {isstring} = "GaAs"
    plot_options.which_plot (1,1) string = "All"
    plot_options.disp_permittivity (1,:) {mustBeInteger} = []
    plot_options.disp_truncation (1,:) {mustBeInteger} = []
    plot_options.plot_layers (1,:) {mustBeInteger} = []
    plot_options.field_direc = []
    plot_options.field_slice = []
    plot_options.field_bars = []
end

for j = 1:numel(folderName)

    % Loads parameter set
    p = load("results/" + folderName(j) + "/param.mat", "param");
    param = p.param;

    % Loop through parameter set
    for i = 1:numel(param)
        A = load("results/" + folderName(j) + "/sim" + string(i) + ".mat");
        layer = fill_layer(param(i), "results/" + folderName(j));

        % find gaas entry (option which layer)
        gaaspos = find([layer.material] == plot_options.which_layer_jsc, 1);
        if isempty(gaaspos)
            warning("Requested material not found, setting requested jsc to first layer")
            jsc(i) = A.fom(1);
        else
            jsc(i) = A.fom(2)+A.fom(3);%A.fom(gaaspos);
        end

        % displayDiscretized(layer(1).geometry.eps_struc, 2)

        if lower(plot_options.which_plot) == "abs"
            pl = 1;
        elseif lower(plot_options.which_plot) == "haze"
            pl = 2;
        elseif lower(plot_options.which_plot) == "harmonics"
            pl = 3;
        elseif lower(plot_options.which_plot) == "all"
            pl = 4;
        else
            error("Wrong input")
        end

        if isempty(plot_options.whichsims) || ismember(i, plot_options.whichsims)
            RCWA_plot(param(i), A.Sz, layer, i, "Sim " + string(i), pl)

            % Plot permittivity
            if ismember(i, plot_options.disp_permittivity)
                [~] = import_permittivities({layer.material}, param(i).wavelengthArray, true, i);
            end

            % Display truncation scheme
            if ismember(i, plot_options.disp_truncation)
                disp_trunc(param(i).Hmax, param(i).tr_ind, i)
            end

            if param(i).calcFields
                show_fields(param, A.fields, plot_options.field_direc, plot_options.field_slice, plot_options.field_bars)
            end

            if ismember(i, plot_options.plot_layers)
                for jj = 1:numel(layer)
                    if iscell(layer(jj).input)
                        plotLayer(layer, jj)
                    end
                end
            end
        end

    end
    
    if numel(param) > 1
    figure
    tiledlayout('flow')
    for i = 1:numel(param(1).bk)
        for k = 1:numel(param(1).bk(i).list)
            jsk_bk(k) = jsc(param(1).bk(i).list(k));
            val(k) = param(param(1).bk(i).list(k)).(param(1).bk(i).name);
        end
        nexttile
        plot(val, jsk_bk)
        title("Jsc per simulation", "FontWeight", 'bold')
        xname = param(1).bk(i).name;
        xname = string([upper(xname{1}(1)), xname{1}(2:end)]);
        xlabel(xname, "FontWeight", 'bold')
        ylabel("Jsc (mA/cm^2)", "FontWeight", 'bold')
    end
    end

end

end

function disp_trunc(Hmax, tr_ind, sim_nr)
M = - (Hmax - 1) / 2:(Hmax - 1) / 2;
N = - (Hmax - 1) / 2:(Hmax - 1) / 2;
TMAP = zeros(Hmax,Hmax);
TMAP(tr_ind) = 1;

figure
imagesc(TMAP);
title("Truncation scheme, sim " + sim_nr)
xticks(1:Hmax)
yticks(1:Hmax)
xticklabels(M)
yticklabels(N)
end


function plotLayer(layer,i)
figure
if layer(i).reverse
    mesh(sum(layer(i).geometry.eps_struc .* reshape(layer(i).L,1,1,[]),3))
    set(gca, 'zdir', 'reverse')
    zt = get(gca, 'ZTick');
    set(gca, 'ZTickLabel', fliplr(zt))
else
    mesh(sum(layer(i).geometry.eps_struc .* reshape(layer(i).L,1,1,[]),3))
end
title("Layer " + i + ": " + layer(i).material)
zlabel("-Z")
end