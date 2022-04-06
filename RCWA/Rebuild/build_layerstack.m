function layer = build_layerstack(layer)

for i=1:numel(layer)

    % If it's not already, change rough layer thickness into an array
    if (layer(i).roughdim>1)&&(numel(layer(i).L)==1)
        layer(i).L = layer(i).L/layer(i).roughdim*ones(1,layer(i).roughdim);
    end

    % Assign uniform layers
    if iscell(layer(i).input)

        for j=1:numel(layer(i).input)
            if isa(layer(i).input{j},'Surface')
                if numel(layer(i).input)==1
                    input(:, :, 1) = layer(i).input.surfMatrix;
                    numLay = numel(layer(i).material);
                    if numLay > 2
                        warning("Layer is a multilayer but only in the first layer the surface is specified. Proceeding with other layers set to constant thickness.")
                        sz = size(layer(i).input.surfMatrix);
                        firstlaymax = max(layer(i).input.surfMatrix, [], 'all');
                        laythick = sum(layer(i).L) - firstlaymax;
                        if laythick < 0
                            warning("First layer already exceeds maximum layer thickness as specified in layers.xlsx; setting constant layers to 0.")
                            laythick = max(laythick, 0);
                        end
                        if numLay == 3
                            laythick = max(sum(layer(i).L) - firstlaymax, 0);
                            input(:, :, 2) = laythick * ones(sz);
                        elseif numLay >= 4
                            laythick = laythick / (numLay - 2);
                            for j = 2:numLay - 1
                                input(:, :, j) = laythick * ones(sz);
                            end
                        end
                    end
                else
                    input(:,:,j) = layer(i).input{j}.surfMatrix;
                end
            else
                error("Invalid layer input.")
            end
        end

        % Discretize
        [Z, Lnew, Lrecalc] = discretize_surface(input, layer(i).roughdim, ...
            layer(i).tolerance, layer(i).reverse, layer(i).optimRough, layer(i).fill, layer(i).add);

        layer(i).geometry.eps_struc = Z;

        if layer(i).recalcRoughL
            layer(i).L = Lrecalc/layer(i).roughdim*ones(1,layer(i).roughdim);
        end

        if layer(i).optimRough
            layer(i).L = Lnew' * sum(layer(i).L)/layer(i).roughdim;
        end

    elseif layer(i).input == 0 && numel(layer(i).material) == 1
        layer(i).geometry.eps_struc = 1;
        layer(i).geometry.mu  = 1;
    elseif layer(i).input == 0 && numel(layer(i).material) > 1
        warning("Multilayer has no surface profile assigned, resetting to uniform.")
        layer(i).geometry.eps_struc = ones(1, 1, numel(layer(i).material));
        layer(i).geometry.mu  = ones(1, 1, numel(layer(i).material));
    else
        error("Invalid layer input.")
    end


end
end

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
numLayers = N3+1;
eps = 1:numLayers;
y = zeros(NX, NY, N3);

% Minimum 0
for i = 1:N3
    y(:, :, i) = input(:, :, i) - min(input(:, :, i), [], 'all');
    if i>1
        layersUpToNow = sum(y(:, :, 1:i-1), 3);
        vq = interpLayerMax(layersUpToNow);
        y(:, :, i) = add(i-1) * input(:, :, i) + fill(i-1) * max(vq - layersUpToNow,0);
    end
end

totmax = max(sum(y, 3), [], 'all');
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
V(V == 0) = eps(j + 1);

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



% function Z_out = create_shape(t,res_x,res_y)
%
% Z_in = zeros(res_x,res_y);
%
% if t==1
%     % Grating
%     pgon = polyshape([0 0 0.5*res_y 0.5*res_y], [0 res_x res_x 0]);
% elseif t==2
%     % Grating
%     pgon = polyshape([0 0 res_y res_y], [0 0.5*res_x 0.5*res_x 0]);
% elseif t==3
%     % Grating
%     pgon = polyshape([0 0 0.5*res_y 0.5*res_y], [0 0.5*res_x 0.5*res_x 0]);
% elseif t==4
%     % Triangle
%     pgon = polyshape([0.1*res_y 0.5*res_y 0.9*res_y], [0.1*res_x 0.9*res_x 0.1*res_x]);
% elseif t==5
%     %Circle
%     n=20;
%     theta = (0:n-1)*(2*pi/n);
%     r = 0.35*res_x;
%     xc = 0.5*res_y;
%     yc = 0.5*res_x;
%     x = xc + r*cos(theta);
%     y = yc + r*sin(theta);
%     pgon = polyshape(x,y);
% end
% [colGrid,rowGrid] = meshgrid(1:size(Z_in,2),1:size(Z_in,1));
% idx = isinterior(pgon,[colGrid(:),rowGrid(:)]);
% idx = reshape(idx,size(Z_in));
%
% Z_out = Z_in;
% Z_out(idx) = 1;
% end

% function Z_out = roughsurf(Xres,Yres,Zres)
%
% rng(1337)
% N = [Xres Yres];
% F = 2;
% [X,Y] = ndgrid(1:N(1),1:N(2));
% i = min(X-1,N(1)-(X-1));
% j = min(Y-1,N(2)-(Y-1));
%
% H = exp(-.5*(i.^2+j.^2)/F^2);
% Z = real(ifft2(H.*fft2(randn(N))));
%
% %normalize
% data=Z/max(max(abs(Z)));
%
% d=discretize(data,linspace(-1,1,Zres));
%
% Z_out=zeros(N(1),N(2),Zres);
% for i=1:Zres
%     Z_out(:,:,i)=i<=d;
% end
%
% end
