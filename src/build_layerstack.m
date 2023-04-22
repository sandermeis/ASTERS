function layer = build_layerstack(layer, param)

% Loop through layers
for i = 1:numel(layer)

    %% Structured layers

    % Check if layer is cell, else check for uniform, else error invalid input
    if iscell(layer(i).input)
        % If L is a scalar, and layer is not uniform, change layer thickness
        % into a vector
        if (layer(i).roughdim > 1) && (numel(layer(i).L) == 1)
            layer(i).L = layer(i).L / layer(i).roughdim * ones(1, layer(i).roughdim);
        end

        % Loop through composite layer
        for j = 1:numel(layer(i).input)
            %% Surface object
            % Check if layer is a Surface object, else error invalid input
            if isa(layer(i).input{j},'Surface') && layer(i).roughdim > 1
                % If layer is a multilayer, but layer cell only has one Surface entry, else fill as normal
                if numel(layer(i).input) == 1
                    input(:, :, 1) = layer(i).input{1}.surfMatrix;
                    numLay = numel(layer(i).material);
                    % if number of materials smaller than number of multilayers
                    if numLay == 1
                        warning("Multilayer has 2 components, however only 1 material(s) are specified, setting to previous layer material")
                        if i == 1
                            if layer(i).reverse
                                layer(i).material(2) = param.ref_medium;
                            else
                                layer(i).material(2) = layer(i).material(1);
                                layer(i).material(1) = param.ref_medium;
                            end
                        else
                            if layer(i).reverse
                                layer(i).material(2) = layer(i - 1).material(end);
                            else
                                layer(i).material(2) = layer(i).material(1);
                                layer(i).material(1) = layer(i - 1).material(end);
                            end
                        end
                    elseif numLay > 2
                        warning("Layer is a multilayer but only in the first layer the surface is specified. Proceeding with other layers set to constant thickness.")
                        sz = size(layer(i).input{1}.surfMatrix);
                        firstlaymax = max(layer(i).input{1}.surfMatrix, [], 'all');
                        % Thickness left in multilayer
                        laythick = sum(layer(i).L) - firstlaymax;
                        if laythick < 0
                            warning("First layer already exceeds maximum layer thickness as specified in layers.xlsx; setting constant layers to 0.")
                            laythick = max(laythick, 0);
                        end
                        % If 2 layers input is unchanged
                        % If 3 layers, add a uniform layer to input with
                        % thickness to bring max to total defined thickness
                        if numLay == 3
                            laythick = max(sum(layer(i).L) - firstlaymax, 0);
                            input(:, :, 2) = laythick * ones(sz);
                        % If 4 or more layers divide remaining thickness evenly over uniform layers
                        elseif numLay >= 4
                            laythick = laythick / (numLay - 2);
                            for k = 2:numLay - 1
                                input(:, :, k) = laythick * ones(sz);
                            end
                        end
                    end
                %% Multiple surface inputs
                else
                    input(:, :, j) = layer(i).input{j}.surfMatrix;
                end

                % Discretize
                [Z, Lnew, Lrecalc] = discretize_surface(input, layer(i).roughdim, ...
                    param.tolerance, layer(i).reverse, param.optimRough, param.fill, param.add);

                layer(i).geometry.eps_struc = Z;

                if param.recalcRoughL
                    layer(i).L = Lrecalc / layer(i).roughdim * ones(1, layer(i).roughdim);
                    warning(sprintf("Recalculated %s layer thickness to %d nm", join(layer(i).material(:), ", "), Lrecalc))
                end

                if param.optimRough
                    layer(i).L = Lnew' * sum(layer(i).L) / layer(i).roughdim;
                end
            %% Already discretized input
            % Check if there is 1 input, and it is numeric, square and binary
            elseif (numel(layer(i).input) == 1) && isnumeric(layer(i).input{1}) && (size(layer(i).input{1}, 1) == size(layer(i).input{1}, 2)) && all((layer(i).input{1} == 0) | (layer(i).input{1} == 1), 'all')
                warning("Assuming binary input")
                layer(i).geometry.eps_struc = layer(i).input{1};
                layer(i).geometry.eps_struc(layer(i).geometry.eps_struc == 1) = 2;
                layer(i).geometry.eps_struc(layer(i).geometry.eps_struc == 0) = 1;

                if layer(i).roughdim ~= size(layer(i).input{j}, 3)
                    layer(i).roughdim = size(layer(i).input{j}, 3);
                    layer(i).L = sum(layer(i).L) / layer(i).roughdim * ones(1, layer(i).roughdim);
                    warning("Input array does not match roughdim, overwriting roughdim")
                end
                
                if numel(layer(i).material) ~= max(layer(i).geometry.eps_struc, [], 'all')
                    % numer of materials smaller than number of multilayers
                    if numel(layer(i).material) < max(layer(i).geometry.eps_struc, [], 'all')
                        warning("Multilayer has 2 components, however only 1 material(s) are specified, setting to previous layer material")
                        if i == 1 % Possibly need to add something for reverse
                            layer(i).material(2) = layer(i).material(1);
                            layer(i).material(1) = param.ref_medium;
                        else
                            layer(i).material(2) = layer(i).material(1);
                            layer(i).material(1) = layer(i - 1).material(end);
                        end
                        % number of materials larger than number of multilayers
                    else
                        error("Too many materials for input multilayer")
                    end
                end
            else
                error("Invalid layer input.")
            end
        end

    %% Uniform layers
    elseif layer(i).input == 0 && numel(layer(i).material) == 1
        layer(i).geometry.eps_struc = 1;
        layer(i).geometry.mu  = 1;
        if layer(i).roughdim > 1
            warning("Roughdim can't be larger than 1 for uniform layers, overwriting roughdim to 1")
            layer(i).roughdim = 1;
        end
    %% Uniform multilayer
    elseif layer(i).input == 0 && numel(layer(i).material) > 1
        warning("Multilayer has no surface profile assigned, resetting to evenly spaced uniform multilayer")
        num_mat = numel(layer(i).material);
        layer(i).geometry.eps_struc = ones(1, 1, num_mat);
        layer(i).geometry.mu  = 1;
        layer(i).L = sum(layer(i).L) / num_mat * ones(1, num_mat);
    else
        error("Invalid layer input.")
    end

end
end


function [dnew, Lnew, Lrecalc] = discretize_surface(Z, Zres, eps, rev, opt, fill, add)

[Z_out, Lrecalc] = discretizeLayers(Z, Zres, fill, add);

if rev
    Z_out = flip(Z_out, 3);
end
% Optimize probably doesnt work with reverse
if opt
    [dnew, Lnew] = optim_discr(Z_out, eps);
else
    dnew = Z_out;
    Lnew = 1;
end

end


function [V, totmax] = discretizeLayers(input, Zres, fill, add)

[NX, NY, N3] = size(input);

% Multilayer numbering
numLayers = N3 + 1;
eps = 1:numLayers;

y = zeros(NX, NY, N3);

% Loop through multilayer surfaces
for i = 1:N3
    % Set minimum to zero
    y(:, :, i) = input(:, :, i) - min(input(:, :, i), [], 'all');
    % If multilayer
    if i>1
        % Get sum of underlying stack of layers
        layersUpToNow = sum(y(:, :, 1:i - 1), 3);

        % Interpolate
        vq = interpLayerMax(layersUpToNow);

        % Add/fill next layer
        y(:, :, i) = add(i - 1) * input(:, :, i) + fill(i - 1) * max(vq - layersUpToNow, 0);
    end
end
% Peak of multilayer
totmax = max(sum(y, 3), [], 'all');

% Z grid
dz = totmax / Zres;
z = linspace(dz, totmax, Zres);

V = zeros(NX, NY, Zres);
f_cond = zeros(NX, NY, numLayers);

for i = 1:Zres
    for j = 1:numLayers - 1
        f_cond(:, :, j + 1) = (z(i) <= sum(y(:, :, 1:j), 3));

        f(:, :, j) = eps(j) * (~f_cond(:, :, j) & f_cond(:, :, j + 1));
    end
    V(:, :, i) = sum(f, 3);
end
% Also assign "rest" surface
V(V == 0) = eps(end);

end


function vq = interpLayerMax(v_in)

[NX, NY] = size(v_in);

% Extend boundaries to account for PBC
superv = [v_in, v_in, v_in;...
        v_in, v_in, v_in;...
        v_in, v_in, v_in];

% Maximum in both rows and columns
supermaxima = islocalmax(superv, 1) & islocalmax(superv, 2);
% Grid with size of superv
[X_in, Y_in] = ndgrid(1:size(supermaxima, 1), 1:size(supermaxima, 2));
% Interpolate just between the maxima
F = scatteredInterpolant(X_in(supermaxima), Y_in(supermaxima), superv(supermaxima));
% Recover old size grid
[Xnew, Ynew] = ndgrid(NX+1:2*NX, NY+1:2*NY);
vq = F(Xnew, Ynew);

end


function [dnew, Lnew] = optim_discr(d, eps)
[Nx, Ny, Nz] = size(d);
sameEdges = zeros(1, Nz);
sameEdges(1) = 1;
k = 1;
d_st = d(:, :, 1);
for i = 2:length(sameEdges)

    likeness = d(:, :, i) == d_st;
    if (1 - sum(likeness(:)) / (Nx * Ny)) > eps %larger than margin, so doesnt match
        k = k + 1;
        d_st = d(:, :, i);
    end
    sameEdges(i) = k;
end

[~, ia, ic] = unique(sameEdges);
dnew = d(:, :, ia);
Lnew = accumarray(ic, 1);

end
