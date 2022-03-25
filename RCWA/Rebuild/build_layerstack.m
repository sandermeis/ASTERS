function layer = build_layerstack(layer)

for i=1:numel(layer)

    % Make an array of layer thickness
    if (layer(i).roughdim>1)&&(numel(layer(i).L)==1)
        layer(i).L = layer(i).L/layer(i).roughdim*ones(1,layer(i).roughdim);
    end

    % uniform layers
    if layer(i).input==0
        layer(i).geometry.eps_struc = 1;
        layer(i).geometry.mu  = 1;
    else
        %optional optimize, optional eps, optional layer height, feedback
        if isa(layer(i).input,'Surface')
            % check for correct dimensions
            input = layer(i).input.surfMatrix;
            [Z, Lnew, Lrecalc] = discretize_surface(input, layer(i).roughdim, ...
                layer(i).tolerance,'reverse',layer(i).reverse,'optimize',layer(i).optimRough);

            layer(i).geometry.eps_struc = Z;

            if layer(i).recalcRoughL
                layer(i).L = Lrecalc/layer(i).roughdim*ones(1,layer(i).roughdim);
            end

            if layer(i).optimRough
                layer(i).L = Lnew' * sum(layer(i).L)/layer(i).roughdim;
            end
        else
            error("Layer not a Surface object")
            %[Z, L, Lrecalc] = realsurf(input,layer(i).roughdim, param.optimRough, param.reverse, param.tolerance);
        end
    end
end
end

function [dnew, Lnew, Lrecalc] = discretize_surface(Z, Zres, eps, options)
arguments
    Z (:,:) {mustBeNumeric}
    Zres uint8 {mustBeNumeric}
    eps (1,1) {mustBeNumeric}
    options.reverse logical = false
    options.optimize logical = false
end

Z = Z - min(Z(:));

[Nx, Ny] = size(Z);

Lrecalc = max(Z(:));

d = discretize(Z,Zres+1)-1;

Z_out = zeros(Nx, Ny, Zres);

for i = 1:Zres
    Z_out(:,:,i) = (d>=i);
end

if options.reverse
    Z_out = flip(Z_out,3);
end
% Optimize probably doesnt work with reverse, should work only for nonzero value,
%currently also counts zeros
if options.optimize
    [dnew,Lnew] = optim_discr(Z_out, eps);
else
    dnew = Z_out;
    Lnew = 1;
end

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
