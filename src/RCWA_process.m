function RCWA_process(folderName, plot_options)
arguments
    folderName (1,:) string
    plot_options.whichsims (1,:) {mustBeInteger} = []
    plot_options.which_layer_jsc (1,1) {isstring} = "GaAs"
    plot_options.which_plot (1,1) string = "All"
    plot_options.plot_4d (1,:) {mustBeInteger} = []
    plot_options.disp_permittivity (1,:) {mustBeInteger} = []
    plot_options.disp_truncation (1,:) {mustBeInteger} = []
    plot_options.plot_surface (1,:) {mustBeInteger} = []
    plot_options.plot_discretized (1,:) {mustBeInteger} = []
    plot_options.avg = []
    plot_options.field_direc {mustBeMember(plot_options.field_direc,["x", "y", "z", "norm", "abs"])} = "norm"
    plot_options.field_slice {mustBeMember(plot_options.field_slice,["x", "y", "z"])} = "x"
    plot_options.field_bars = true;
end

jsc_plot = [];

for j = 1:numel(folderName)

    % Loads parameter set
    p = load("results/" + folderName(j) + "/param.mat", "param");
    param = p.param;

    % Loop through parameter set
    for i = 1:numel(param)

        A = load("results/" + folderName(j) + "/sim" + string(i) + ".mat");
        layer = fill_layer(param(i), "results/" + folderName(j));
        Sz{i} = A.Sz;
        % Calculate short circuit current
        fom = jsc(squeeze(sum(A.Sz, 1)), param(i).wavelengthArray);

        % find gaas entry (option which layer)
        if numel(param(i).wavelengthArray) > 1
            gaaspos = find([layer.material] == plot_options.which_layer_jsc, 1);
            if isempty(gaaspos)
                warning("Requested material not found, setting requested jsc to first layer")
                jsc_plot(i) = fom(1);
            else
                jsc_plot(i) = fom(gaaspos);
            end
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
            RCWA_plot(param(i), A.Sz, layer, "Sim " + string(i), pl)

            % Plot permittivity
            if ismember(i, plot_options.disp_permittivity)
                [~] = import_permittivities({layer.material}, param(i).wavelengthArray, true, i);
            end

            % Display truncation scheme
            if ismember(i, plot_options.disp_truncation)
                disp_trunc(param(i).Hmax, param(i).tr_ind, i)
            end

            % Plot of fields
            if param(i).calcFields
                show_fields(param, A.fields, plot_options.field_direc, plot_options.field_slice, plot_options.field_bars)
            end

            % Plot non discretized surface
            if ismember(i, plot_options.plot_surface)
                for jj = 1:numel(layer)
                    if iscell(layer(jj).input)
                        plotLayer(layer, jj)
                    end
                end
            end

            % Plot discretized surface
            if ismember(i, plot_options.plot_discretized)
                for jj = 1:numel(layer)
                    if iscell(layer(jj).input)
                        displayDiscretized(layer(jj).geometry.eps_struc)
                    end
                end
            end

            % Plot 4d harmonics
            if ismember(i, plot_options.plot_4d)
                tile_harmonics4d(param, Sz{i}, layer, "Sim " + string(i))
            end
        end

    end

    if numel(param) > 1
        figure
        tiledlayout('flow')
        for i = 1:numel(param(1).bk)
            for k = 1:numel(param(1).bk(i).list)
                jsk_bk(k) = jsc_plot(param(1).bk(i).list(k));
                val(k) = param(1).bk(i).val(k);
            end
            nexttile
            plot(val, jsk_bk)
            title("Jsc per simulation", "FontWeight", 'bold')
            xname = param(1).bk(i).name;
            xname = string([upper(xname{1}(1)), xname{1}(2:end)]);
            xlabel(xname, "FontWeight", 'bold')
            ylabel("Jsc (mA/cm^2)", "FontWeight", 'bold')
        end

        if ~isempty(plot_options.avg)
            avg_list_nr = find(plot_options.avg==[param(1).bk.name]);
            avg_sz_cell = Sz(param(1).bk(avg_list_nr).list);
            avg_sz = avg_sz_cell{1};
            for items = 2:numel(avg_sz_cell)
                avg_sz = avg_sz + avg_sz_cell{items};
            end
            avg_sz = avg_sz / numel(avg_sz_cell);

            RCWA_plot(param(param(1).bk(avg_list_nr).list(1)), avg_sz, layer, "Average over parameter " + plot_options.avg, pl)
        end
    end

end

end

function n = fixLayerString(n)
% Add R and T to layer string
n(end+1) = {"R"};
n(end+1) = {"T"};
for j=1:numel(n)
    n{j} = [n{j}{:}];
end
end


function disp_trunc(Hmax, tr_ind, sim_nr)
M = - (Hmax - 1) / 2:(Hmax - 1) / 2;
N = - (Hmax - 1) / 2:(Hmax - 1) / 2;
TMAP = zeros(Hmax, Hmax);
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
    mesh(sum(layer(i).geometry.eps_struc .* reshape(layer(i).L, 1, 1, []), 3))
    set(gca, 'zdir', 'reverse')
    zt = get(gca, 'ZTick');
    set(gca, 'ZTickLabel', fliplr(zt))
else
    mesh(sum(layer(i).geometry.eps_struc .* reshape(layer(i).L, 1, 1, []), 3))
end
title("Layer " + i + ": " + layer(i).material)
zlabel("-Z")
end


function RCWA_plot(param, Sz, layer, titlestring, whichdisp)
switch whichdisp
    case 1
        plot_results(param, Sz, layer, titlestring);
    case 2
        plothaze(param, Sz, layer, titlestring);
    case 3
        jsc_harmonics(param, Sz, layer, titlestring);
    case 4
        plot_results(param, Sz, layer, titlestring);
        plothaze(param, Sz, layer, titlestring);
        jsc_harmonics(param, Sz, layer, titlestring);
end
end


function show_fields(param, fields, direc, slice, bars)

direc_nr = find(direc == ["x", "y", "z", "norm", "abs"]);

layer = fill_layer(param);
boundaries = [layer.L];
wl = 1;
%wavelength % material % component
if direc_nr == 4
    d = real(abs(fields{wl}{1}) + abs(fields{wl}{2}) + abs(fields{wl}{3}));
    title_text = "|E|";
elseif direc_nr == 5
    d = real(abs(fields{wl}{1}) + abs(fields{wl}{2}) + abs(fields{wl}{3}));
    title_text = "P_{abs}";
else
    d = real(abs(fields{wl}{direc_nr}));
    title_text = "|E" + direc + "|";
end


if direc_nr ~= 5

    plot_field_3d(d, param, layer, slice, bars, boundaries, title_text)

else
    c = 1;%2.998e8;%2*0.002654418728*256;
    for iter = 1:numel(param.wavelengthArray)

        for i = 1:numel(layer)

            layer(i).geometry.eps = layer(i).geometry.eps_struc;
            for j = 1:numel(layer(i).permittivities)
                % Assign permittivity at specific wavelength
                eps = layer(i).permittivities{j}(param.wavelengthArray(iter)); % wl using interpolant
                layer(i).geometry.eps(layer(i).geometry.eps_struc == j) = eps;
            end
        end
        k = 1;
        for i = 1:numel(layer)
            for j = 1:numel(layer(i).L)
                for l = 1:layer(i).L(j)
                    absorp(:, :, k) = -0.5 * 2 * pi / param.wavelengthArray(iter) * d(:, :, k) .* imag(layer(i).geometry.eps(:, :, j));
                    k = k + 1;
                end
            end
        end

    end
    plot_field_3d(absorp, param, layer, slice, bars, boundaries, title_text)

end
end


function plot_field_3d(d, param, layer, slice, bars, boundaries, title_text)
cmax = max(d, [], 'all');
cmin = min(d, [], 'all');
switch slice
    case {"x", "y"}
        figure
        for ii = 1:size(d, 2)
            imagesc(squeeze(d(:, ii, :)))
            colorbar
            set(gca, 'clim', [cmin cmax])
            set(gca, 'YDir','normal')
            set(gca, 'XTick', [0:0.1:1] * size(d, 3), 'XTickLabel', [0:0.1:1] * sum(boundaries)) % 10 ticks
            set(gca, 'YTick', [0:0.1 * param.res:param.res], 'YTickLabel', [0:0.1 * param.size / param.res:param.size / param.res] * param.res) % 20 ticks
            xlabel("z (nm)")
            ylabel(slice + " (nm)")
            title(title_text + ", " + slice + " = " + (ii - 1) + " nm")
            for iii = 1:numel(layer)
                h{iii} = text(sum([layer(1:iii).L]) - sum(layer(iii).L) / 2, 0.01 * param.res, layer(iii).material);
                set(h{iii}, 'Rotation', 90);
            end

            if bars
                for i = 1:numel(boundaries)
                    xline(sum(boundaries(1:i)))
                end
            end
            drawnow
            pause(0.1)
        end

    case "z"

        figure
        for ii = 1:size(d, 3)
            imagesc(squeeze(d(:, :, ii)))
            colorbar
            set(gca, 'clim', [cmin cmax])
            set(gca, 'YDir','normal')
            set(gca, 'XTick', [0:0.1 * param.res:param.res], 'XTickLabel', [0:0.1 * param.size / param.res:param.size / param.res] * param.res)
            set(gca, 'YTick', [0:0.1 * param.res:param.res], 'YTickLabel', [0:0.1 * param.size / param.res:param.size / param.res] * param.res) % 20 ticks
            xlabel("x (nm)")
            ylabel("y (nm)")
            title(title_text + ", z = " + (ii-1) + " nm")

            drawnow
            pause(0.1)
        end
end
end


function plothaze(param, Sz, layer, titlestring)

centralH = squeeze(abs(Sz((end + 1) / 2, :, :)));
sumH = squeeze(sum(abs(Sz), 1));
diffH = sumH - centralH;
haze = arrayfun(@(a, b) a / b * (b > 1e-12), diffH, sumH);

n = fixLayerString({layer.material});

figure
colororder(parula(size(haze, 1)))
for i = 1:size(haze, 1)
    hold on
    plot(param.wavelengthArray, abs(haze(i, :)).', 'LineWidth', 2)
end
xlabel('Wavelength (nm)')
ylabel('Haze')
legend(n, 'location', 'eastoutside')
title(titlestring, "FontSize", 16, "FontWeight", 'bold')
end


function plot_results(param, Sz_in, layer, titlestring)

plot_type = "area";
Sz = squeeze(sum(Sz_in, 1)).';
N = size(Sz, 2);
x_grid = repmat(param.wavelengthArray', 1, N);
n = fixLayerString({layer.material});

% Place GaAs on bottom
gaaspos = find(n == "GaAs" | n == "H/GaAs" | n == "GaAs_1meV" | n == "GaAs_5meV"...
    | n == "GaAs_10meV" | n == "GaAs_15meV" | n == "GaAs_20meV" | n == "GaAs_25meV", 1);
if ~isempty(gaaspos)
    Sz(:, [1 gaaspos]) = Sz(:, [gaaspos 1]);
    n([1 gaaspos]) = n([gaaspos 1]);
end

% Calculate jsc
for i = 1:N
    jsc_plot(i) = jsc(Sz(:, i)', param.wavelengthArray);
end

h = figure('Color', 'w', 'Position', [100 100 1300 600]);

colors = [76, 144, 186; 43, 194, 194; 244, 184, 17; 222, 102, 62; 255, 145, 43] / 255;

colororder(colors);
hx = subplot(1, 4, [1, 2, 3]);

if plot_type=="lines"
    ar = plot(x_grid, Sz, 'LineWidth', 2);
    title(titlestring, "FontSize", 18, "FontWeight", 'bold')
    xlim([param.wavelengthArray(1), param.wavelengthArray(end)])
    ylim([0, 1])
    xlabel("Wavelength (nm)", "FontSize", 16, "FontWeight", 'bold')
    ylabel("Absorption (a.u.)", "FontSize", 16, "FontWeight", 'bold')
    set(hx,'FontSize', 14, 'LineWidth', 2)
else
    ar = area(x_grid, Sz, 'LineWidth', 2);
    title(titlestring, "FontSize", 18, "FontWeight", 'bold')
    xlim([param.wavelengthArray(1), param.wavelengthArray(end)])
    ylim([0, 1])
    xlabel("Wavelength (nm)", "FontSize", 16, "FontWeight", 'bold')
    ylabel("Absorption (a.u.)", "FontSize", 16, "FontWeight", 'bold')
    set(hx, 'FontSize', 14, 'LineWidth', 2)

    % Jsc box
    legstring = pad(n) + sprintf(" Jsc: ") + pad(compose('%0.2f', string(jsc_plot)), 'left') + sprintf(" mA/cm^2");
    hx2 = subplot(1, 4, 4);
    hx2.Visible = 'off';
    [~,legend_h,~,~] = legendflex(ar, cellstr(legstring), 'ref', hx2, 'anchor', [1 1], 'buffer', [0 0], 'FontName', 'monospaced');


    % Hatched area for many layers

    hstyle = {'single', 'single', 'cross', 'single', 'single', 'cross'};
    hdir = {0, 45, 0, 90, 135, 45};

    for i = 1:numel(ar)
        h_ind = ceil(i / length(colors)) - 1;
        if h_ind > 0 && h_ind < length(colors)
            j = mod(i - 1, length(hstyle)) + 1;
            hatchfill2(ar(i), hstyle{j}, 'HatchAngle', hdir{j}, 'HatchDensity', 40, 'HatchColor', 'k', 'HatchLineWidth', 2)
            hatchfill2(legend_h(length(ar) + i), hstyle{j}, 'HatchAngle', hdir{j}, 'HatchDensity', 40, 'HatchColor', 'k', 'HatchLineWidth', 2)
        end
    end
end

end


function jsc_harmonics(param, Sz, layer, titlestring)
% Plot jsc contribution for every harmonic
n = fixLayerString({layer.material});

num_lay = size(Sz, 2);
h = figure;
t = tiledlayout(h, 'flow');
title(t, "Jsc per harmonic per material, " + titlestring, "FontSize", 16, "FontWeight", 'bold')
cmin_old = 0;
cmax_old = 0;
for i = 1:num_lay
    jsc_empty = zeros(param.P, param.Q);
    h(i) = nexttile;
    jsc_empty(param.tr_ind) = jsc(squeeze(Sz(:, i, :)), param.wavelengthArray);
    im = imagesc(jsc_empty);
    title(n(i))
    set(h(i), "Color", "none", 'YDir', 'normal', 'TickLength', [0 0])
    im.AlphaData = jsc_empty ~= 0;
    xticks(1:param.P)
    yticks(1:param.Q)
    xticklabels(h(i), num2cell(-0.5 * (param.P - 1):0.5 * (param.P - 1)))
    yticklabels(h(i), num2cell(-0.5 * (param.Q - 1):0.5 * (param.Q - 1)))

    h2(i) = axes(t);
    h2(i).Layout.Tile = h(i).Layout.Tile;
    set(h2(i), "Color", "none", 'YDir', 'normal')
    h2(i).XLim = [0, param.P];
    h2(i).YLim = [0, param.Q];
    h2(i).YTickLabel = {};
    h2(i).XTickLabel = {};
    xticks(h2(i), 0:param.P)
    yticks(h2(i), 0:param.Q)
    h2(i).XGrid = 'on';
    h2(i).YGrid = 'on';
    h2(i).GridAlpha = 1;
    j_nocentral = jsc_empty;
    j_nocentral(round(param.P / 2), round(param.Q / 2)) = NaN;
    cmin_new = min(j_nocentral, [], 'all');
    cmax_new = max(j_nocentral, [], 'all');
    cmin_old = max(cmin_old, cmin_new);
    cmax_old = max(cmax_old, cmax_new);
end
set(h, 'Colormap', jet, 'CLim', [cmin_old cmax_old])
cbh = colorbar(h(end));
cbh.Layout.Tile = 'east';

end


function tile_harmonics4d(param, Sz, layer, titlestring)
% Display all harmonics for every wavelength
P = param.P;
Q = param.Q;
tr_ind = param.tr_ind;

n = fixLayerString({layer.material});

h = figure;
t = tiledlayout(h, 'flow', 'TileSpacing', 'Compact');
title(t, titlestring)
cmax = max(Sz(:));
cmin = 0;
a = repair_harmonics(Sz, P, Q, tr_ind);

[X, Y, Z] = ndgrid(-((P - 1) / 2):((P - 1) / 2),-((Q - 1) / 2):((Q - 1) / 2), param.wavelengthArray);

for i = 1:size(Sz, 2)
    D = squeeze(a(:, :, i, :));
    nexttile
    scatter3(Z(:), X(:), Y(:), 400, D(:), 'filled');
    xlabel("Wavelength (nm)")
    ylabel("P")
    zlabel("Q")
    title(n(i))
    axis ij
    caxis manual
    caxis([cmin cmax]);
    alpha color
    alpha scaled
end
cb = colorbar;
cb.Layout.Tile = 'east';

end


function a = repair_harmonics(Sz, P, Q, tr_ind)
% Set truncated harmonics back to untruncated
num_lay = size(Sz, 2);
num_lab = size(Sz, 3);
a = zeros(P, Q, num_lay, num_lab);
for i = 1:num_lay
    for j = 1:num_lab
        c = a(:, :, i, j);
        c(tr_ind) = Sz(:, i, j);
        a(:, :, i, j) = c;
    end
end
end


function displayDiscretized(V)
% Displays discretized layer in blocks
numLayers = max(V, [], 'all');
clrmap = hsv(numLayers);
for i = 1:numLayers
    figure
    for j = 1:i
        patch(FindExternalVoxels(V == j), 'FaceColor', clrmap(j,:), 'LineWidth', 0.1, 'FaceAlpha', 1/j, 'EdgeAlpha', 0.5)
        hold on
    end
    view(45, 30)
    grid on
end
end