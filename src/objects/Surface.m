classdef Surface < handle
    %SURFACE Produces a surface, to which features can be added.
    %   s = SURFACE(r), if r is a scalar, produces a Surface with resolution r
    %   and default size.
    %
    %   s = SURFACE(r, sz) specifies the size sz of the Surface in nanometers.
    %
    %   Surface properties:
	%              surfMatrix - SURFACE grid matrix.
    %                Features - Feature struct list.
    %                surfsize - Size of the SURFACE in nanometers.
    %                 surfres - Resolution of the SURFACE in pixels.
    %                   surfx - x coordinates of the surface vector.
    %                   surfy - y coordinates of the surface vector.
    %                   surfX - x coordinates of the surface matrix.
    %                   surfY - y coordinates of the surface matrix.
    %
    %	Surface methods:
    %              addFeature - Adds Feature object to the SURFACE.
    %            listFeatures - List of Features that have been added.
    %                    plot - Displays current SURFACE grid.
    %            clearSurface - Clears current SURFACE grid.
    %            placeFeature - Places specified Feature on the grid.
    %                  rotate - Rotates SURFACE grid.
    %       addRandomFeatures - Adds specified Feature to random locations.
    %                  hscale - Scales current SURFACE grid.
    %            addRoughsurf - Adds specified rough surface to SURFACE grid.
    %           placeFeatures - Place selected Features on the SURFACE grid.
    
    %   Copyright 2022 Sander Meis.
    
    properties
        %SURFMATRIX - Magnitude of the surface matrix.
        %   surfMatrix is a matrix giving the magnitude of the SURFACE.
        %
        %   See also SURFACE
        surfMatrix
        
        %FEATURES - List of features that have been added to the SURFACE.
        %   Features is a list of structs .
        %
        %   See also SURFACE
        Features (1, :) struct
    end
    properties (SetAccess = protected)
        %SURFSIZE - Size of the SURFACE in nanometers.
        %   Currently Surfaces can only be square, so SURFSIZE represents both
        %   the x and y length. 
        %   On construction, the size can be chosen explicitly, else it
        %   defaults to 10000.
        %
        %   See also SURFACE
        surfsize = 10000
        
        %SURFRES - Resolution of the SURFACE in pixels.
        %   Currently Surfaces can only be square, so SURFRES represents both
        %   the x and y resolution. 
        %   On construction, the resolution can be chosen explicitly, else it
        %   defaults to 512.
        %
        %   See also SURFACE
        surfres = 512
        
        %SURFX - x coordinates vector of the Surface.
        %   SURFX is a vector giving the x coordinates of the SURFACE.
        %
        %   See also SURFACE
        surfx
        
        %SURFY - y coordinates vector of the Surface.
        %   SURFY is a vector giving the y coordinates of the SURFACE.
        %
        %   See also SURFACE
        surfy
        
        %SURFX - x coordinates matrix of the Surface.
        %   SURFX is a matrix giving the x coordinates of the SURFACE.
        %
        %   See also SURFACE
        surfX
        
        %SURFY - y coordinates matrix of the Surface.
        %   SURFY is a matrix giving the y coordinates of the SURFACE.
        %
        %   See also SURFACE
        surfY
    end
    
    methods
        function obj = Surface(varargin)
            if nargin == 1
                obj.surfres = varargin{1};
            elseif nargin == 2
                obj.surfres = varargin{1};
                obj.surfsize = varargin{2};
            else
                error("Too many input arguments")
            end
            obj.surfx = linspace(0, obj.surfsize, obj.surfres);
            obj.surfy = linspace(0, obj.surfsize, obj.surfres);
            [obj.surfX, obj.surfY] = ndgrid(obj.surfx, obj.surfy);
            obj.surfMatrix = zeros(obj.surfres, obj.surfres);
        end
        
        
        function obj = addFeature(obj, fobj, xpos, ypos, PBC)
            arguments
                obj
                fobj Feature
                xpos (1, :) {mustBeInteger}
                ypos (1, :) {mustBeInteger}
                PBC logical = true
            end
            
            %ADDFEATURE Adds Feature object to the Surface.
            %   Input is ADDFEATURE(Feature, xpos, ypos, PBC)
            %
            %   Feature - Feature object.
            %   xpos    - X grid position(s) of feature(s)
            %   ypos    - Y grid position(s) of feature(s)
            %   PBC     - Periodic boundary conditions (0 or 1)
            %
            %   See also: ADDFEATURE, PLACEFEATURE, ADDRANDOMFEATURES
            for i = 1:numel(xpos)
                if numel(fobj)==numel(xpos)
                    fobj_index = i;
                else
                    fobj_index = 1;
                end
                if (fobj(fobj_index).resolution <= obj.surfres)
                    if PBC || checkPlacement(obj, xpos(i), ypos(i), fobj(fobj_index).resolution, fobj(fobj_index).resolution)
                        if isempty(obj.Features)
                            obj.Features(1).fobj = fobj(fobj_index);
                            obj.Features(1).xpos = xpos(i);
                            obj.Features(1).ypos = ypos(i);
                        else
                            obj.Features(end + 1).fobj = fobj(fobj_index);
                            obj.Features(end).xpos = xpos(i);
                            obj.Features(end).ypos = ypos(i);
                        end
                    else
                        warning("Position is out of bounds, no object added")
                    end
                else
                    fobj(fobj_index).Z = fobj(fobj_index).Z(1:obj.surfres, 1:obj.surfres);
                    fob(fobj_index).resolution = obj.surfres;
                    warning("Feature larger than surface, Feature sampled down to surface size")
                    if PBC || checkPlacement(obj, xpos(i), ypos(i), fobj(fobj_index).resolution, fobj(fobj_index).resolution)
                        if isempty(obj.Features)
                            obj.Features(1).fobj = fobj(fobj_index);
                            obj.Features(1).xpos = xpos(i);
                            obj.Features(1).ypos = ypos(i);
                        else
                            obj.Features(end + 1).fobj = fobj(fobj_index);
                            obj.Features(end).xpos = xpos(i);
                            obj.Features(end).ypos = ypos(i);
                        end
                    else
                        warning("Position is out of bounds, no object added")
                    end
                end
            end
        end
        
        
        function listFeatures(obj)
            %LISTFEATURES Prints a list with all Features.
            %   Function has no input.
            %
            %   See also: REPORT, ADDFEATURE
            if ~isempty(obj.Features)
                fprintf("Surface ""%s"" consists of these Features:\n", inputname(1))
                fprintf("--------------------\n")
                for i = 1:numel(obj.Features)
                    %number
                    % coordinates
                    % type
                    nr = length(obj.Features(i).xpos);
                    fprintf("--------------------\n")
                    fprintf("Feature %u\nType: ""%s""; Number: %u; Number density: %d num/um^2\nCoordinates: \n", i, obj.Features(i).fobj.shape, nr, nr / (obj.surfsize / 1000)^2)
                    fprintf("Copy %u: (%u, %u)\n", [1:length(obj.Features(i).xpos); obj.Features(i).xpos; obj.Features(i).ypos])
                    fprintf("--------------------\n")
                end
            else
                warning("Surface has no features")
            end
        end
        
        
        function report(obj)
            %REPORT Reports on some surface parameters.
            %   Function has no input. Prints the Variance, RMS roughness,
            %   Skewness, Kurtosis
            %
            %
            %   See also: CLEARFEATURES, CLEARSURFACE, ADDFEATURE
            fprintf("Surface ""%s"" Report:\n", inputname(1))
            fprintf("--------------------\n")
            meansurf = mean(obj.surfMatrix(:));
            %Mean roughness
            fprintf("meansurf: %d\n", meansurf)
            var = mean((obj.surfMatrix(:) - meansurf).^2);
            fprintf("var: %d\n", var)
            rms = sqrt(mean((obj.surfMatrix(:) - meansurf).^2));
            fprintf("rms: %d\n", rms)
            skew = mean((obj.surfMatrix(:) - meansurf).^3 ./ rms.^3);
            fprintf("skew: %d\n", skew)
            kurt = mean((obj.surfMatrix(:) - meansurf).^4 ./ rms.^4);
            fprintf("kurt: %d\n", kurt)
            
            figure
            tiledlayout('flow')
            %Normalized height distribution
            nexttile
            histogram(obj.surfMatrix(:) / numel(obj.surfMatrix), 'EdgeColor', 'none');
            title("Height distribution")

            %Power spectral density
            nexttile
            m = obj.surfres;
            a = obj.surfsize / obj.surfres;
            qx_1 = zeros(m, 1);
            for k = 0:m - 1
                qx_1(k + 1) = (2 * pi / m) * (k);
            end
            qx_2 = fftshift(qx_1);
            qx_3 = unwrap(qx_2 - 2 * pi);
            qx = qx_3 / a;
            surf(qx, qx, fft2(obj.surfMatrix) .* conj(fft2(obj.surfMatrix)), 'EdgeColor', 'none');
            title("PSD")

            %Radially averaged power spectral density
            nexttile
            [q , C, ~] = psd_2D(obj.surfMatrix, obj.surfsize / obj.surfres);
            loglog(q, C)
            title("Radially averaged PSD")
            shading flat

            %Autocorrelation function
            nexttile
            surf(abs(fftshift(ifft2(fft2(obj.surfMatrix) .* conj(fft2(obj.surfMatrix))))) ./ (obj.surfres.^2))
            title("Autocorrelation function")
            shading interp
        end
        
        
        function placeFeature(obj, i, xpos, ypos, PBC, mode)
            arguments
                obj
                i
                xpos
                ypos
                PBC logical = true
                mode string = "add"
            end
            %PLACEFEATURE Places all Features on the grid.
            %   Input is PLACEFEATURE(i, xpos, ypos, options)
            %
            %   i       - Feature list number.
            %   xpos    - X grid position(s) of feature(s)
            %   ypos    - Y grid position(s) of feature(s)
            %
            %   Where options can be the following string:
            %   'PBC'   - Periodic boundary conditions (0 or 1).
            %   'mode'  - "add", "merge" or "replace".
            %
            %   See also: ADDFEATURE, PLACEFEATURE, ADDRANDOMFEATURES
            
            nx = obj.Features(i).fobj.resolution;
            ny = obj.Features(i).fobj.resolution;
            for pos = 1:numel(xpos)
                x = xpos(pos);
                y = ypos(pos);
                if PBC
                    ind = ((y: y + ny - 1) .* ones(ny, 1) - 1) * obj.surfres * 2 + (x: x + nx - 1)' .* ones(1, nx);
                    k = Surface.toPBC(obj.surfres * 2, obj.surfres, obj.surfres, ind(:));
                    temp = zeros(obj.surfres);
                    temp(k) = obj.Features(i).fobj.Z;
                    if mode == "add"
                        obj.surfMatrix = obj.surfMatrix + temp;
                    elseif mode == "merge"
                        temp2 = temp - obj.surfMatrix;
                        obj.surfMatrix = obj.surfMatrix + (temp2 + abs(temp2)) / 2;
                    elseif mode == "replace"
                        obj.surfMatrix(k) = obj.Features(i).fobj.Z;
                    else        %exclude mode?, only place when no other features in the way
                        error('wrong mode')
                    end
                else
                    if checkPlacement(obj, x, y, nx, ny)
                        obj.surfMatrix(x:x + nx - 1, y:y + ny - 1) = obj.surfMatrix(x: x + nx - 1, y: y + ny - 1) + obj.Features(i).fobj.Z;
                    else
                        errorString = 'Placement of Feature ' + string(i) + ...
                            ' at position x = ' + string(x) + ' of ' + string(obj.surfres) + ...
                            ', y = ' + string(y) + ' of ' + string(obj.surfres) + ', with size ΔX = ' + ...
                            string(nx) + ', ΔY = ' + string(ny) + ' failed, skipping.';
                        warning(errorString)
                    end
                end
            end
        end
        
        
        function placeFeatures(obj, options)
            arguments
                obj
                options.PBC logical = true
                options.mode string = "add"
            end
            %PLACEFEATURES Places all Features on the grid.
            %   Input is PLACEFEATURES(options)
            %
            %   Where options can be the following string:
            %   'PBC'    	- Periodic boundary conditions (0 or 1).
            %   'mode'      - "add", "merge" or "replace".
            %
            %   See also: ADDFEATURE, PLACEFEATURE, ADDRANDOMFEATURES
            for i = 1:numel(obj.Features)
                x = obj.Features(i).xpos;
                y = obj.Features(i).ypos;
                placeFeature(obj, i, x, y, options.PBC, options.mode)
            end
        end
        
        
        function addRandomFeatures(obj, fobj, number, options)
            arguments
                obj
                fobj
                number
                options.PBC logical = true
                options.seed = 0
                %options.mode string = "add"
            end
            %ADDRANDOMFEATURES Adds a Feature object with random coordinates to the Surface a certain number of times.
            %   Input is ADDRANDOMFEATURES(Feature, number, options)
            %
            %   Feature     - Feature object.
            %   number      - How many features are to be added.
            %
            %   Where option can be the following string:
            %   'PBC'    	- Periodic boundary conditions (0 or 1).
            %
            %   See also: ADDFEATURE, PLACEFEATURE, ADDRANDOMFEATURES
            rng(options.seed)
            if ~options.PBC
                nx = fobj.resolution;
                ny = fobj.resolution;
                maxx = size(obj.surfMatrix, 1) - nx + 1;
                maxy = size(obj.surfMatrix, 2) - ny + 1;
                
                if numel(fobj.height) > 1
                    for i = 1:numel(fobj.height)
                        fcopy = fobj;
                        fcopy.height = fobj.height(i);
                        fcopy.Z = fcopy.height * fcopy.Z / max(fcopy.Z, [], 'all');
                        xpos = randi([1, maxx], 1, 1);
                        ypos = randi([1, maxy], 1, 1);
                        addFeature(obj, fcopy, xpos, ypos, options.PBC);
                    end
                else
                    xpos = randi([1, obj.surfres], 1, number);
                    ypos = randi([1, obj.surfres], 1, number);
                    addFeature(obj, fobj, xpos, ypos, options.PBC);
                end
            else
                if numel(fobj.height) > 1
                    for i = 1:numel(fobj.height)
                        fcopy = fobj;
                        fcopy.height = fobj.height(i);
                        fcopy.Z = fcopy.height * fcopy.Z / max(fcopy.Z, [], 'all');
                        xpos = randi([1, obj.surfres], 1, 1);
                        ypos = randi([1, obj.surfres], 1, 1);
                        addFeature(obj, fcopy, xpos, ypos, options.PBC);
                    end
                else
                    xpos = randi([1, obj.surfres], 1, number);
                    ypos = randi([1, obj.surfres], 1, number);
                    addFeature(obj, fobj, xpos, ypos, options.PBC);
                end
            end
        end
        
        %%% Placeholder
        function addRandomFeaturesMC(obj, fobj, tarRMS, number, options)
            arguments
                obj
                fobj
                tarRMS
                number
                options.PBC logical = true
                %options.mode string = "add"
            end
            if ~options.PBC
                nx = fobj.resolution;
                ny = fobj.resolution;
                maxx = size(obj.surfMatrix, 1) - nx + 1;
                maxy = size(obj.surfMatrix, 2) - ny + 1;
                
                xpos = randi([1, maxx], 1, number);
                ypos = randi([1, maxy], 1, number);
                %placeFeature(obj, i, x, y, options.PBC, options.mode);
                addFeature(obj, fobj, xpos, ypos, options.PBC);
                
            else
                for i = 1:number
                xpos = randi([1, obj.surfres], 1, number);
                ypos = randi([1, obj.surfres], 1, number);
                %placeFeature(obj, i, x, y, options.PBC, options.mode);
                surfmean = mean(obj.surfMatrix(:));
                rms = sqrt(mean((obj.surfMatrix(:) - surfmean).^2));
                delta_rms = tarRMS - rms;
                if delta_rms > 0
                    addFeature(obj, fobj, xpos, ypos, options.PBC);
                end
                end
            end
        end
        
        
        function clearSurface(obj)
            %CLEARSURFACE Clears the surface, this does not clear the Features.
            %   Function has no input
            %
            %   See also: CLEARFEATURES, ADDFEATURE, PLACEFEATURE
            obj.surfMatrix = zeros(obj.surfres, obj.surfres);
        end
        
        
        function clearFeatures(obj)
            %CLEARFEATURES Clears all Features, this does not clear the surface matrix.
            %   Function has no input
            %
            %   See also: CLEARSURFACE, ADDFEATURE, PLACEFEATURE
            obj.Features = struct([]);
        end
        
        
        function clearAll(obj)
            %CLEARALL Clears the surface and Features.
            %   Function has no input
            %
            %   See also: CLEARFEATURES, CLEARSURFACE, ADDFEATURE
            clearSurface(obj)
            clearFeatures(obj)
        end
        
        
        function resize(obj, surfres_new, surfsize_new)
            %RESIZE Resamples Surface to new resolution and size.
            %   Input is RESIZE(res, size)
            %
            %   res     - New resolution.
            %   size    - New size in nm.
            %
            %   See also: HSCALE, PLACEFEATURE, ADDRANDOMFEATURES
            clearFeatures(obj)
            surfx_new = linspace(0, surfsize_new, surfres_new);
            surfy_new = linspace(0, surfsize_new, surfres_new);
            [surfX_new, surfY_new] = ndgrid(surfx_new, surfy_new);
            F = griddedInterpolant({obj.surfx, obj.surfy}, obj.surfMatrix);
            obj.surfMatrix = F({surfx_new, surfy_new});
            obj.surfres = surfres_new;
            obj.surfsize = surfsize_new;
            obj.surfx = surfx_new;
            obj.surfy = surfy_new;
            obj.surfX = surfX_new;
            obj.surfY = surfY_new;
        end
        
        
        function hscale(obj, h)
            %HSCALE Rescales the Surface.
            %   Input is HSCALE(h)
            %
            %   h       - New height.
            %
            %   See also: ADDFEATURE, PLACEFEATURE, ADDRANDOMFEATURES
            omax = max(obj.surfMatrix(:));
            if omax > 0
                obj.surfMatrix = h * obj.surfMatrix / omax;
            else
                warning("Can't rescale empty surface")
            end
        end
        
        
        function plot(obj)
            %PLOT Displays Surface object.
            %   Function has no input
            %   See also: ADDFEATURE, PLACEFEATURE, ADDRANDOMFEATURES
            figure
            mesh(obj.surfX, obj.surfY, obj.surfMatrix)
            title("Surface Architecture")
            xlabel("X (nm)")
            ylabel("Y (nm)")
            zlabel("Height (nm)")
        end
        
        
        function obj = addRoughsurf(obj, options)
            arguments
                obj
                options.PBC logical = true
                options.mode string = "Rough"
                options.sigma (1, 1) double = 100
                options.hurst (1, 1) double = 0.5
                options.height (1, 1) double = 100
                options.x (1, 1) {mustBeInteger} = 1
                options.y (1, 1) {mustBeInteger} = 1
                options.seed (1, 1) = 1
            end
            %ADDROUGHSURF Adds a rough layer to the Surface as a Feature object.
            %   Input is ADDROUGHSURF(options)
            %
            %   Where options can be the following strings:
            %   'PBC'    	- Periodic boundary conditions (0 or 1).
            %   'mode'      - "Rough" or "Artificial".
            %   'sigma'  	- Standard deviation or RMS roughness (only for "Artificial").
            %   'hurst'  	- Fractal parameter.
            %   'height'    - Height (only for "Rough").
            %   'x'         - Feature x coordinate.
            %   'y'         - Feature y coordinate.
            %
            %   See also: ADDFEATURE, PLACEFEATURE, ADDRANDOMFEATURES
            if options.mode == "Artificial"
                [Z, ~, ~] = artificial_surf(options.sigma, options.hurst, obj.surfsize, obj.surfres, obj.surfres);
            elseif options.mode == "Rough"
                Z = Surface.roughsurf(obj.surfres, obj.surfres, options.height, 1 / options.hurst, options.seed);
            else
                error("Invalid mode selected")
            end
            addFeature(obj, Feature(Z, obj.surfsize), options.x, options.y, options.PBC);
        end


        function obj = addUniform(obj, t)

            %ADDROUGHSURF Adds a uniform layer to the Surface as a Feature object.
            %   Input is ADDUNIFORM(t)
            %
            %   't'    	    - Thickness.
            %
            %   See also: ADDROUGHSURF
            addFeature(obj, Feature(t * ones(obj.surfres, obj.surfres), obj.surfsize), 1, 1);
        end

        
        function obj = rotate(obj)
            %ROTATE Rotates Surface object 90 degrees anticlockwise.
            %   Function has no input
            %   See also: ADDFEATURE, PLACEFEATURE, ADDRANDOMFEATURES
            obj.surfMatrix = rot90(obj.surfMatrix);
        end
        
    end
    methods (Access = private)
        function allowedPlacement = checkPlacement(obj, x, y, nx, ny)
            if all(round(x) == x) && all(round(y) == y) && all(x >= 1) && all(y >= 1) && all(x + nx - 1 <= obj.surfres) && all(y + ny - 1 <= obj.surfres)
                allowedPlacement = true;
            else
                allowedPlacement = false;
            end
        end
    end
    methods (Static)
        function Z = roughsurf(Xres, Yres, height, F, seed)
            rng(seed)
            N = [Xres Yres];
            [X, Y] = ndgrid(1:N(1), 1:N(2));
            i = min(X - 1, N(1) - (X - 1));
            j = min(Y - 1, N(2) - (Y - 1));
            
            H = exp(-.5 * (i.^2 + j.^2) / F^2);
            Z = real(ifft2(H .* fft2(randn(N))));
            
            Z = Z - min(Z(:));
            Z = height * Z / max(Z(:));
        end


        function lind_pbc = toPBC(NX, NXnew, NYnew, ind)
            r = mod(ind - 1, NX) + 1;
            c = ceil(ind ./ NX);
            r_pbc = mod(r - 1, NXnew) + 1;
            c_pbc = mod(c - 1, NYnew) + 1;
            lind_pbc = (c_pbc - 1) .* NXnew + r_pbc;
        end

    end
end

function mustBeEqualSize(a, b)
if ~isequal(size(a), size(b))
    eid = 'Size:notEqual';
    msg = 'Size of inputs must be equal.';
    throwAsCaller(MException(eid, msg))
end
end