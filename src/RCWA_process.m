function batch = RCWA_process(folderName, plot_options)
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

batch = struct;
jsc_plot = [];
mymap = readmatrix("src/data/cm_coldwarm.csv");

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

        % Plot absorption, haze and harmonics
        if isempty(plot_options.whichsims) || ismember(i, plot_options.whichsims)
            plot_out = RCWA_plot(param(i), A.Sz, layer, "Sim " + string(i), pl, mymap);

             f = fieldnames(plot_out);
             for fi = 1:length(f)
                batch(j).data(i).(f{fi}) = plot_out.(f{fi});
             end
             
        end

        % Plot permittivity
        if ismember(i, plot_options.disp_permittivity)
            [~] = import_permittivities({layer.material}, param(i).wavelengthArray, true, i, mymap);
        end

        % Display truncation scheme
        if ismember(i, plot_options.disp_truncation)
            disp_trunc(param(i).Hmax, param(i).tr_ind, i)
        end

        % Plot of fields
        if param(i).calcFields
            show_fields(param, A.fields, plot_options.field_direc, plot_options.field_slice, plot_options.field_bars)
            batch(j).data(i).fields = A.fields;
        end

        % Plot surface as surf
        if ismember(i, plot_options.plot_surface)
            for jj = 1:numel(layer)
                if iscell(layer(jj).input)
                    plotLayer(param, layer, jj, mymap)
                end
            end
        end

        % Plot surface as blocks
        if ismember(i, plot_options.plot_discretized)
            for jj = 1:numel(layer)
                if iscell(layer(jj).input)
                    displayDiscretized(param, layer, jj, mymap)
                end
            end
        end

        % Plot 4d harmonics
        if ismember(i, plot_options.plot_4d)
            tile_harmonics4d(param, Sz{i}, layer, "Absorption per harmonic per wavelength, Sim " + string(i), mymap)
        end

    % Return data
    batch(j).data(i).param = param(i);
    batch(j).data(i).layer = layer;
    batch(j).data(i).Sz = A.Sz;
    batch(j).data(i).jsc = fom;
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
            plot(val, jsk_bk, "Color", "k", LineWidth=2)
            xname = param(1).bk(i).name;
            xname = string([upper(xname{1}(1)), xname{1}(2:end)]);
            title("Jsc per " + xname, "FontSize", 16, "FontWeight", "normal", "Interpreter", "none")
            xlabel(xname, "FontSize", 14, "Interpreter", "none")
            ylabel("Jsc (mA/cm^2)", "FontSize", 14)

            batch(j).("series_" + xname).x = val;
            batch(j).("series_" + xname).y = jsk_bk;
        end

        if ~isempty(plot_options.avg)
            avg_list_nr = find(plot_options.avg==[param(1).bk.name]);
            if ~isempty(avg_list_nr)
                avg_sz_cell = Sz(param(1).bk(avg_list_nr).list);
                avg_sz = avg_sz_cell{1};
                for items = 2:numel(avg_sz_cell)
                    avg_sz = avg_sz + avg_sz_cell{items};
                end
                avg_sz = avg_sz / numel(avg_sz_cell);
    
                batch(j).("avg_" + plot_options.avg) = RCWA_plot(param(param(1).bk(avg_list_nr).list(1)), avg_sz, layer, "averaged over parameter " + plot_options.avg, pl, mymap);
                batch(j).("avg_" + plot_options.avg).Sz = avg_sz;
            else
                error("Unable to find parameter " + plot_options.avg)
            end
        end
    end

end

end


function n = fixLayerString(n)
% Add R and T to layer string
n(end+1) = {"R"};
n(end+1) = {"T"};
for j = 1:numel(n)
    n{j} = strjoin(n{j}, ", ");
end
n = cellstr(n);
end


function disp_trunc(Hmax, tr_ind, sim_nr)
M = - (Hmax - 1) / 2:(Hmax - 1) / 2;
N = - (Hmax - 1) / 2:(Hmax - 1) / 2;
TMAP = zeros(Hmax, Hmax);
TMAP(tr_ind) = 1;

figure
colormap(gray)
imagesc(~TMAP);
title("Truncation scheme, sim " + sim_nr, "FontSize", 16, "FontWeight", "normal")
xticks(1:Hmax)
yticks(1:Hmax)
xticklabels(M)
yticklabels(N)
xlabel("M", "FontSize", 14)
ylabel("N", "FontSize", 14)
end


function plotLayer(param, layer, i, mymap)
figure
grid = linspace(0, param.size, param.res);
colormap(mymap)
[X, Y] = ndgrid(grid);
if layer(i).reverse
    surf(X, Y, sum((layer(i).geometry.eps_struc == 1) .* reshape(layer(i).L, 1, 1, []), 3), 'FaceColor', 'interp', 'EdgeColor', 'interp')
    set(gca, 'zdir', 'reverse')
    zt = get(gca, 'ZTick');
    set(gca, 'ZTickLabel', fliplr(zt))
else
    surf(X, Y, sum((layer(i).geometry.eps_struc == 1) .* reshape(layer(i).L, 1, 1, []), 3), 'FaceColor', 'interp', 'EdgeColor', 'interp')
end
title("Layer " + i + ": " + layer(i).material(1), "FontSize", 16, "FontWeight", "normal")
xlim('tight'); ylim('tight'); zlim('tight')
colorbar
axis equal
xlabel("X", "FontSize", 14)
ylabel("Y", "FontSize", 14)
zlabel("Z", "FontSize", 14)
end


function displayDiscretized(param, layer, jj, mymap)
V = layer(jj).geometry.eps_struc;
% Displays discretized layer in blocks
numLayers = max(V, [], 'all');
clrmap = mymap(round(linspace(1, 1024, numLayers)), :);
for i = 1:numLayers
    f = figure;
    title("Layer " + jj + ": " + strjoin("" + layer(jj).material, ", "), "FontSize", 16, "FontWeight", "normal")
    xlabel("X", "FontSize", 14)
    ylabel("Y", "FontSize", 14)
    zlabel("Z", "FontSize", 14)
    for j = 1:i
        patch(FindExternalVoxels(V == j), 'FaceColor', clrmap(j,:), 'LineWidth', 0.1, 'FaceAlpha', 1/j, 'EdgeAlpha', 0.5)
        hold on
    end
    if layer(jj).reverse
        f.CurrentAxes.ZDir = 'Reverse';
    end
    xlim('tight'); ylim('tight'); zlim('tight')
    axis equal
    view(45, 30)
    grid on
end
end


function data = RCWA_plot(param, Sz, layer, titlestring, whichdisp, mymap)

switch whichdisp
    case 1
        data.absorption = plot_results(param, Sz, layer, "Absorption per layer, " + titlestring, mymap);
    case 2
        data.haze = plothaze(param, Sz, layer, "Haze per layer, " + titlestring, mymap);
    case 3
        data.jsc_harm = jsc_harmonics(param, Sz, layer, "Jsc per harmonic per layer, " + titlestring, mymap);
    case 4
        data.absorption = plot_results(param, Sz, layer, "Absorption per layer, " + titlestring, mymap);
        data.haze = plothaze(param, Sz, layer, "Haze per material, " + titlestring, mymap);
        data.jsc_harm = jsc_harmonics(param, Sz, layer, "Jsc per harmonic per layer, " + titlestring, mymap);
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
            xlabel("z (nm)", "FontSize", 14)
            ylabel(slice + " (nm)", "FontSize", 14)
            title(title_text + ", " + slice + " = " + (ii - 1) + " nm", "FontSize", 16, "FontWeight", "normal")
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
            xlabel("x (nm)", "FontSize", 14)
            ylabel("y (nm)", "FontSize", 14)
            title(title_text + ", z = " + (ii-1) + " nm", "FontSize", 16, "FontWeight", "normal")

            drawnow
            pause(0.1)
        end
end
end


function haze = plothaze(param, Sz, layer, titlestring, mymap)

centralH = squeeze(abs(Sz((end + 1) / 2, :, :)));
sumH = squeeze(sum(abs(Sz), 1));
diffH = sumH - centralH;
haze = arrayfun(@(a, b) a / b * (b > 1e-12), diffH, sumH);

n = fixLayerString({layer.material});

figure
colororder(mymap(round(linspace(1, 1024, size(haze, 1))), :))
for i = 1:size(haze, 1)
    hold on
    plot(param.wavelengthArray, abs(haze(i, :)).', 'LineWidth', 2)
end
xlabel('Wavelength (nm)', "FontSize", 14)
ylabel('Haze', "FontSize", 14)
legend(n, 'location', 'eastoutside')
title(titlestring, "FontSize", 16, "FontWeight", "normal", "Interpreter", "none")
end


function Sz = plot_results(param, Sz_in, layer, titlestring, mymap)

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

%colors = [76, 144, 186;...
%          43, 194, 194;...
%          244, 184, 17;...
%          222, 102, 62;...
%          255, 145, 43] / 255;

colors = mymap([128, 320, 512, 704, 896], :);

colororder(colors);
hx = subplot(1, 4, [1, 2, 3]);

if plot_type=="lines"
    ar = plot(x_grid, Sz, 'LineWidth', 2);
    title(titlestring, "FontSize", 16, "FontWeight", "normal", "Interpreter", "none")
    xlim([param.wavelengthArray(1), param.wavelengthArray(end)])
    ylim([0, 1])
    xlabel("Wavelength (nm)", "FontSize", 14)
    ylabel("Absorption (a.u.)", "FontSize", 14)
    set(hx,'FontSize', 14, 'LineWidth', 2)
else
    ar = area(x_grid, Sz, 'LineWidth', 2);
    title(titlestring, "FontSize", 16, "FontWeight", "normal", "Interpreter", "none")
    xlim([param.wavelengthArray(1), param.wavelengthArray(end)])
    ylim([0, 1])
    xlabel("Wavelength (nm)", "FontSize", 14)
    ylabel("Absorption (a.u.)", "FontSize", 14)
    set(hx, 'FontSize', 14, 'LineWidth', 2)

    % Jsc box
    legstring = pad(n) + sprintf(" Jsc: ") + pad(compose('%0.2f', string(jsc_plot)), 'left') + sprintf(" mA/cm^2");
    hx2 = subplot(1, 4, 4);
    hx2.Visible = 'off';
    [~,legend_h,~,~] = legendflex(ar, cellstr(legstring), 'ref', hx2, 'anchor', [1 1], 'buffer', [-55 -5], 'FontName', 'monospaced');


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


function output = jsc_harmonics(param, Sz, layer, titlestring, mymap)
% Plot jsc contribution for every harmonic
n = fixLayerString({layer.material});

num_lay = size(Sz, 2);
h = figure;
t = tiledlayout(h, 'flow');
title(t, titlestring, "FontSize", 16, "FontWeight", "normal", "Interpreter", "none")
cmin_old = 0;
cmax_old = 0;
for i = 1:num_lay
    jsc_empty = zeros(param.P, param.Q);
    h(i) = nexttile;
    jsc_empty(param.tr_ind) = jsc(squeeze(Sz(:, i, :)), param.wavelengthArray);
    im = imagesc(jsc_empty);
    output{i} = jsc_empty;
    title(n(i), "FontSize", 14, "FontWeight", "normal")
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
set(h, 'Colormap', mymap, 'CLim', [cmin_old cmax_old])
cbh = colorbar(h(end));
cbh.Layout.Tile = 'east';

end


function tile_harmonics4d(param, Sz, layer, titlestring, mymap)
% Display all harmonics for every wavelength
P = param.P;
Q = param.Q;
tr_ind = param.tr_ind;

n = fixLayerString({layer.material});

h = figure;
t = tiledlayout(h, 'flow', 'TileSpacing', 'Compact');
title(t, titlestring, "FontSize", 16, "FontWeight", "normal")
colormap(mymap)
a = repair_harmonics(Sz, P, Q, tr_ind);
a(round(P/2), round(Q/2), :, :) = 0;
a = abs(a);
a(a < 1 / (P * Q)) = 0;
a = a / max(a(:));

[X, Y, Z] = ndgrid(-((P - 1) / 2):((P - 1) / 2),-((Q - 1) / 2):((Q - 1) / 2), param.wavelengthArray);

for i = 1:size(Sz, 2)
    D = squeeze(a(:, :, i, :));
    nexttile
    scatter3(Z(:), X(:), Y(:), 200, D(:), 'filled');
    xlabel("Wavelength (nm)", "FontSize", 14)
    ylabel("P", "FontSize", 14)
    zlabel("Q", "FontSize", 14)
    title(n(i), "FontSize", 14, "FontWeight", "normal")
    axis ij
    caxis manual
    caxis([0 1]);
    alpha color
    alpha scaled
    alpha('none')
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