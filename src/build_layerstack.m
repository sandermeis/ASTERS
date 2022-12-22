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
                    if numLay > 2
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
                    % something off, third index probs doesnt work with
                    % discretize
                    input(:,:,j) = layer(i).input{j}.surfMatrix;
                end

                % Discretize
                [Z, Lnew, Lrecalc] = discretize_surface(input, layer(i).roughdim, ...
                    param.tolerance, layer(i).reverse, param.optimRough, param.fill, param.add);

                layer(i).geometry.eps_struc = Z;

                if param.recalcRoughL
                    layer(i).L = Lrecalc / layer(i).roughdim * ones(1, layer(i).roughdim);
                    warning(sprintf("Recalculated %s layer thickness to %d nm", join(layer(i).material(:), ", "), Lrecalc))
                end

                if layer(i).optimRough
                    layer(i).L = Lnew' * sum(layer(i).L) / layer(i).roughdim;
                end
            %% Already discretized input
            % Check if there is 1 input, and it is numeric, square and binary
            elseif (numel(layer(i).input) == 1) && isnumeric(layer(i).input{1}) && (size(layer(i).input{1},1)==size(layer(i).input{1},2)) && all((layer(i).input{1}==0)|(layer(i).input{1}==1),'all')
                warning("Assuming binary input")
                layer(i).geometry.eps_struc = layer(i).input{1};
                layer(i).geometry.eps_struc(layer(i).geometry.eps_struc==1)=2;
                layer(i).geometry.eps_struc(layer(i).geometry.eps_struc==0)=1;

                if layer(i).roughdim ~= size(layer(i).input{j},3)
                    layer(i).roughdim = size(layer(i).input{j},3);
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
        layer(i).geometry.mu  = ones(1, 1, num_mat);
        layer(i).L = sum(layer(i).L) / num_mat * ones(1, num_mat);
    else
        error("Invalid layer input.")
    end

end
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
