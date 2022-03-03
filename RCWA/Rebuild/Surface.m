classdef Surface < handle
    
    properties
        surfMatrix
        Features (1,:) struct
    end
    properties (SetAccess = immutable)
        surfsize = 10000
        surfres = 512
        surfx
        surfy
        surfX
        surfY
    end
    
    methods
        function obj = Surface(varargin)
            if nargin==1
                obj.surfsize = varargin{1};
            elseif nargin>=2
                obj.surfsize = varargin{1};
                obj.surfres = varargin{2};
            end
            obj.surfx = linspace(0, obj.surfsize, obj.surfres);
            obj.surfy = linspace(0, obj.surfsize, obj.surfres);
            [obj.surfX, obj.surfY] = meshgrid(obj.surfx, obj.surfy);
            obj.surfMatrix = zeros(obj.surfres, obj.surfres);
        end
        
        
        function obj = addFeature(obj,fobj,xpos,ypos,PBC)
            arguments
                obj
                fobj Feature
                xpos (1,:) {mustBeInteger}
                ypos (1,:) {mustBeInteger}
                PBC logical = true
            end
            if PBC||checkPlacement(obj,xpos,ypos,fobj.resolution,fobj.resolution)
                if isempty(obj.Features)
                    obj.Features(1).fobj = fobj;
                    obj.Features(1).xpos = xpos;
                    obj.Features(1).ypos = ypos;
                else
                    obj.Features(end+1).fobj = fobj;
                    obj.Features(end).xpos = xpos;
                    obj.Features(end).ypos = ypos;
                end
            else
                warning("Position is out of bounds, no object added")
            end
            
        end
        
        
        function listFeatures(obj)
            if ~isempty(obj.Features)
                fprintf("Surface ""%s"" consists of these Features:\n",inputname(1))
                fprintf("--------------------\n")
                for i=1:numel(obj.Features)
                    %number
                    % coordinates
                    % type
                    nr = length(obj.Features(i).xpos);
                    fprintf("--------------------\n")
                    fprintf("Feature %u\nType: ""%s""; Number: %u; Number density: %d num/um^2\nCoordinates: \n", i, obj.Features(i).fobj.shape, nr, nr/(obj.surfsize/1000)^2)
                    fprintf("Copy %u: (%u, %u)\n",[1:length(obj.Features(i).xpos); obj.Features(i).xpos; obj.Features(i).ypos])
                    fprintf("--------------------\n")
                end
            else
                fprintf("Surface has no features")
            end
        end
        
        
        function placeFeature(obj,i,x,y,PBC,mode)
            arguments
                obj
                i
                x
                y
                PBC logical = true
                mode string = "add" %add, replace, merge
            end
            
            nx = obj.Features(i).fobj.resolution;
            ny = obj.Features(i).fobj.resolution;
            
            if PBC
                ind = ((y: y + ny-1).*ones(ny,1)-1)*obj.surfres*2+(x: x + nx-1)'.*ones(1,nx);
                k = toPBC(obj.surfres*2,obj.surfres,obj.surfres,ind(:));
                temp = zeros(obj.surfres);
                temp(k) = obj.Features(i).fobj.Z;
                if mode=="add"
                    obj.surfMatrix = obj.surfMatrix + temp;
                elseif mode=="merge"
                    temp2 = temp - obj.surfMatrix;
                    obj.surfMatrix = obj.surfMatrix + (temp2+abs(temp2))/2;
                elseif mode=="replace"
                    obj.surfMatrix(k) = obj.Features(i).fobj.Z;
                else        %exclude mode?, only place when no other features in the way
                    error('wrong mode')
                end
            else
                if checkPlacement(obj,x,y,nx,ny)
                    obj.surfMatrix(x: x + nx-1,y: y + ny-1) = obj.surfMatrix(x: x + nx-1,y: y + ny-1) + obj.Features(i).fobj.Z;
                else
                    errorString = 'Placement of Feature ' + string(i) + ...
                        ' at position x = ' + string(x) + ' of ' + string(obj.surfres) + ...
                        ', y = ' + string(y) + ' of ' + string(obj.surfres) + ', with size ΔX = ' + ...
                        string(nx) + ', ΔY = ' + string(ny) + ' failed, skipping.';
                    warning(errorString)
                end
            end
        end
        
        
        function placeFeatures(obj,options)
            arguments
                obj
                options.PBC logical = true
                options.mode string = "add"
            end
            for i=1:numel(obj.Features)
                x = obj.Features(i).xpos;
                y = obj.Features(i).ypos;
                placeFeature(obj,i,x,y,options.PBC,options.mode)
            end
        end
        
        
        function addRandomFeatures(obj,fobj,number,options)
            arguments
                obj
                fobj
                number
                options.PBC logical = true
                %options.mode string = "add"
            end
            if ~options.PBC
                nx = fobj.resolution;
                ny = fobj.resolution;
                maxx = size(obj.surfMatrix,1)-nx+1;
                maxy = size(obj.surfMatrix,2)-ny+1;
                
                xpos = randi([1,maxx],1,number);
                ypos = randi([1,maxy],1,number);
                %placeFeature(obj,i,x,y,options.PBC,options.mode)
                addFeature(obj,fobj,xpos,ypos,options.PBC)
                
            else
                
                xpos = randi([1,obj.surfres],1,number);
                ypos = randi([1,obj.surfres],1,number);
                %placeFeature(obj,i,x,y,options.PBC,options.mode)
                addFeature(obj,fobj,xpos,ypos,options.PBC)
            end
        end
        
        
        function clearSurface(obj)
            obj.surfMatrix = zeros(obj.surfsize*obj.surfres,obj.surfsize*obj.surfres);
        end
        
        
        function hscale(obj,h)
            obj.surfMatrix = obj.surfMatrix*h/max(obj.surfMatrix(:));
        end
        
        
        function plot(obj)
            figure
            mesh(obj.surfX, obj.surfY, obj.surfMatrix)
        end
        
        
        function obj = addRoughsurf(obj,options)
            arguments
                obj
                options.PBC logical = true
                options.mode string = "Rough"
                options.sigma double = 100
                options.hurst double = 0.5
                options.height double = 100
                options.x (1,1) {mustBeInteger} = 1
                options.y (1,1) {mustBeInteger} = 1
            end
            
            if options.mode == "Artificial"
                [Z , ~ , ~] = artificial_surf(options.sigma, options.hurst, obj.surfsize, obj.surfres, obj.surfres);
            elseif options.mode == "Rough"
                Z = Surface.roughsurf(obj.surfres,obj.surfres,options.height,1/options.hurst);
            else
                error("Invalid mode selected")
            end
            addFeature(obj,Feature(Z, obj.surfsize),options.x,options.y,options.PBC);
        end
        
        
        function obj = rotate(obj)
            obj.surfMatrix = rot90(obj.surfMatrix);
        end
        
    end
    methods (Access = private)
        function allowedPlacement = checkPlacement(obj,x,y,nx,ny)
            if all(round(x)==x)&&all(round(y)==y)&&all(x>=1)&&all(y>=1)&&all(x+nx-1<=obj.surfres)&&all(y+ny-1<=obj.surfres)
                allowedPlacement = true;
            else
                allowedPlacement = false;
            end
        end
    end
       methods (Static)
        function Z = roughsurf(Xres,Yres,height,F)
            N = [Xres Yres];
            [X,Y] = ndgrid(1:N(1),1:N(2));
            i = min(X-1,N(1)-(X-1));
            j = min(Y-1,N(2)-(Y-1));
            
            H = exp(-.5*(i.^2+j.^2)/F^2);
            Z = real(ifft2(H.*fft2(randn(N))));
            
            Z = Z-min(Z(:));
            Z = height * Z/max(Z(:));
        end
        
    end
end

