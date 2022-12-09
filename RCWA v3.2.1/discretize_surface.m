function [dnew, Lnew, Lrecalc] = discretize_surface(Z, Zres, eps, rev, opt, fill, add)

%fill = zeros(size(Z, 3)-1);
%add = ones(size(Z, 3)-1);

[Z_out, Lrecalc] = discretizeLayers(Z, Zres, fill, add);

if rev
    Z_out = flip(Z_out, 3);
end
% Optimize probably doesnt work with reverse, should work only for nonzero value,
%currently also counts zeros
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
numLayers = N3+1;
eps = 1:numLayers;

y = zeros(NX, NY, N3);

% Minimum 0
% Loop through multilayer surfaces
for i = 1:N3
    % Set minimum to zero
    y(:, :, i) = input(:, :, i) - min(input(:, :, i), [], 'all');
    % If multilayer
    if i>1
        % Get sum of underlying stack of layers
        layersUpToNow = sum(y(:, :, 1:i-1), 3);

        % Interpolate
        vq = interpLayerMax(layersUpToNow);

        % Add/fill next layer
        y(:, :, i) = add(i-1) * input(:, :, i) + fill(i-1) * max(vq - layersUpToNow,0);
    end
end
% Peak of multilayer
totmax = max(sum(y, 3), [], 'all');

% Z grid
dz = totmax/Zres;
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

superv = [v_in, v_in, v_in;...
    v_in, v_in, v_in;...
    v_in, v_in, v_in];
supermaxima = islocalmax(superv, 1) & islocalmax(superv, 2);
[X_in, Y_in] = ndgrid(1:size(supermaxima, 1), 1:size(supermaxima, 2));
F = scatteredInterpolant(X_in(supermaxima), Y_in(supermaxima), superv(supermaxima));
[Xnew, Ynew] = ndgrid(NX+1:2*NX, NY+1:2*NY);
vq = F(Xnew, Ynew);

end


function [dnew, Lnew] = optim_discr(d, eps)
[Nx, Ny, Nz] = size(d);
sameEdges = zeros(1,Nz);
sameEdges(1) = 1;
k = 1;
d_st = d(:,:,1);
for i = 2:length(sameEdges)

    likeness = d(:,:,i)==d_st;
    if (1-sum(likeness(:))/(Nx*Ny))>eps %larger than margin, so doesnt match
        k = k + 1;
        d_st = d(:,:,i);
    end
    sameEdges(i) = k;
end

[~, ia, ic] = unique(sameEdges);
dnew = d(:,:,ia);
Lnew = accumarray(ic,1);

end

