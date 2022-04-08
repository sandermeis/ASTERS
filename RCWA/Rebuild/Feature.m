classdef Feature
    %FEATURE Produces a specified feature geometry
    %   f = FEATURE(r), if r is a scalar, produces a feature with resolution r
    %   and default size and height.
    %
    %   f = FEATURE(A), if A is a matrix, produces a feature with the
    %   matrix as feature magnitude and default size and height. The input 
    %   matrix must be square.
    %
    %   f = FEATURE("/PATHNAME.csv"), if "/PATHNAME.csv" is a string, produces
    %   a feature with imported data where the specified string is used as
    %   a pathname. The imported matrix must be square.
    %
    %   f = FEATURE(..., sz) specifies the size of the feature in
    %   nanometers. If sz is a vector, create an array of Features with
    %   size corresponding to elements of sz.
    %
    %   f = FEATURE(..., sz, h) specifies the height of the feature in
    %   nanometers. If h is a vector, create an array of Features with
    %   heights corresponding to elements of h.
    %
    %   f = FEATURE(r, sz, h, "Shape") specifies the shape of the feature.
    %   Only possible for non imported data.
    %
    %   Shape can be one of these strings:
    %               "GratingX" - 2D Grating X
    %               "GratingY" - 2D Grating Y
    %              "GratingXY" - 2D Grating XY
    %               "Triangle" - 2D Triangle, 3D Prism
    %                 "Circle" - 2D Circle, 3D Cylinder
    %                 "Sphere" - 3D Half sphere
    %                "Pyramid" - 3D Pyramid
    %                 "RidgeX" - 3D Ridge X
    %                 "RidgeY" - 3D Ridge Y
    %                 "WedgeX" - 3D Wedge X
    %                 "WedgeY" - 3D Wedge Y
    %                   "Cone" - 3D Cone
    %
    %   Feature properties:
    %                    size - Size of the FEATURE in nanometers.
    %              resolution - Resolution of the FEATURE in pixels.
    %                  height - Height of the FEATURE in nanometers.
    %                   shape - Shape of the feature.
    %                       X - x coordinates of the feature matrix.   
    %                       Y - y coordinates of the feature matrix.
    %                       Z - Magnitude of the feature matrix.
    %
    %   Feature methods:
    %                  rotate - Rotates the FEATURE matrix 90 degrees anticlockwise.
    %                    plot - Displays the FEATURE.
    %
    
    %   Copyright 2022 Sander Meis.
    
    properties
        %SIZE - Size of the FEATURE in nanometers.
        %   Currently features can only be square, so SIZE represents both
        %   the x and y length. 
        %   On construction, the size can be chosen explicitly, else it
        %   defaults to 625.
        %
        %   See also FEATURE
        size (1,1) {mustBeNumeric} = 625
        
        %RESOLUTION - Resolution of the FEATURE in pixels.
        %   Currently features can only be square, so RESOLUTION represents both
        %   the x and y resolution. 
        %   On construction, the resolution can be chosen explicitly, else it
        %   defaults to 32.
        %
        %   See also FEATURE
        resolution = 32 %pts/nm
        
        %HEIGHT - Height of the FEATURE in nanometers.
        %   HEIGHT represents the maximum height of the FEATURE.
        %   On construction, the height can be chosen explicitly, else it
        %   defaults to 16.
        %
        %   See also FEATURE
        height (1,:) {mustBeNumeric} = 16
        
        %SHAPE - Shape of the feature.
        %   SHAPE is a string describing the shape of the feature. SHAPE can
        %   be one of:
        %     "GratingX", "GratingY", "GratingXY", "Triangle", "Circle", "Sphere",
        %     "Pyramid", "RidgeX", "RidgeY", "WedgeX", "WedgeY", "Cone", "Custom"
        %
        %   The shape is chosen on construction and cannot be modified. On
        %   construction, the shape can be chosen explicitly, else it
        %   defaults to "Pyramid".
        %
        %   See also FEATURE
        shape string = "Pyramid"
        
    end
    properties (GetAccess = public, SetAccess = ?Surface)
        %X - x coordinates of the feature matrix.
        %   X is a matrix giving the x coordinates of the FEATURE.
        %
        %   See also FEATURE
        X
        
        %Y - y coordinates of the feature matrix.
        %   Y is a matrix giving the y coordinates of the FEATURE.
        %
        %   See also FEATURE
        Y
        
        %Z - Magnitude of the feature matrix.
        %   Z is a matrix giving the magnitude of the FEATURE.
        %
        %   See also FEATURE
        Z
    end
    
    methods
        function obj = Feature(varargin)
            if nargin==0
                warning('No shape added, producing feature with default parameters')
            elseif nargin>=1
                % If single number
                if isnumeric(varargin{1})&&isscalar(varargin{1})
                    obj.resolution = varargin{1};
                    switch nargin
                        case 2
                            obj.size = varargin{2};
                        case 3
                            obj.size = varargin{2};
                            obj.height = varargin{3};
                        case 4
                            obj.size = varargin{2};
                            obj.height = varargin{3};
                            obj.shape = varargin{4};
                        otherwise
                            error('wrong number of input arguments')
                    end
                % If input is a path or a square matrix
                elseif isstring(varargin{1})||ischar(varargin{1})||(isnumeric(varargin{1})&&ismatrix(varargin{1})&&(size(varargin{1},1)==size(varargin{1},2)))
                    
                    if isstring(varargin{1})||ischar(varargin{1})
                        Z = readmatrix(varargin{1});
                    else
                        Z = varargin{1};
                    end
                    N = size(Z);
                    if isnumeric(Z)&&(numel(N)==2)&&(N(1)==N(2))
                        % normalize
                        obj.shape = "Custom";
                        obj.Z = Z;
                        if min(obj.Z, [], 'all') ~= max(obj.Z, [], 'all')
                            obj.Z = obj.Z - min(obj.Z(:));
                        end
                        obj.resolution  = N(1);
 
                        if nargin==2
                            obj.size = varargin{2};

                        % Rescale input height
                        elseif nargin==3
                            obj.size = varargin{2};
                            obj.height = varargin{3};
                            if numel(obj.height)==1
                                obj.Z = obj.height * obj.Z/max(obj.Z(:));
                                warning('Rescaling feature height')
                            end
                        else
                            error('Too many input arguments for custom feature')
                        end
                    else
                        error('Wrong input feature dimensions or data type')
                    end
                else
                    warning('Incorrect input, producing feature with default parameters')
                end
            end
            
            x = linspace(0, obj.size, obj.resolution);
            y = linspace(0, obj.size, obj.resolution);
            [obj.X,obj.Y] = meshgrid(x,y);
            
            if numel(obj.height)>1
                h = 1;
            else
                h = obj.height;
            end

            switch obj.shape
                case "GratingX"
                    % 1D Grating X
                    obj.Z = h * (obj.X<(obj.size/2));
                case "GratingY"
                    % 1D Grating Y
                    obj.Z =  h * (obj.Y<(obj.size/2));
                case "GratingXY"
                    % 2D Grating XY
                    obj.Z =  h * ((obj.X<(obj.size/2))&(obj.Y<(obj.size/2)));
                case "Triangle"
                    % 2D Triangle
                    obj.Z =  h * ((2*abs(obj.X-(obj.size/2))+obj.Y)<obj.size);
                case "Circle"
                    % 2D Circle, 3D Cylinder
                    obj.Z =  h * ((obj.X-obj.size/2).^2+(obj.Y-obj.size/2).^2<=(obj.size/2).^2);
                case "Sphere"
                    % 3D Half sphere
                    tmp = (1-((obj.X-obj.size/2).^2+(obj.Y-obj.size/2).^2)/(obj.size/2)^2);
                    obj.Z =  h * (tmp+abs(tmp))/2;
                case "Pyramid"
                    % 3D Pyramid
                    blaat = @(a,b) max([abs(a-obj.size/2),abs(b-obj.size/2)]);
                    obj.Z =  h * (1-2*arrayfun(blaat,obj.X,obj.Y)/obj.size);
                case "RidgeX"
                    % 3D Ridge X
                    blaat = @(a) max(abs(a-obj.size/2));
                    obj.Z =  h * (1-2*arrayfun(blaat,obj.X)/obj.size);
                case "RidgeY"
                    % 3D Ridge Y
                    blaat = @(a) max(abs(a-obj.size/2));
                    obj.Z =  h * (1-2*arrayfun(blaat,obj.Y)/obj.size);
                case "WedgeX"
                    % 3D Wedge
                    obj.Z =  h * obj.X/max(obj.X(:));
                case "WedgeY"
                    % 3D Wedge
                    obj.Z =  h * obj.Y/max(obj.Y(:));
                case "Cone"
                    tmp = (1-sqrt((obj.X-obj.size/2).^2+(obj.Y-obj.size/2).^2)/(obj.size/2));
                    obj.Z =  h * (tmp+abs(tmp))/2;
                case "Custom"
                    
                otherwise
                    error('Feature name not recognised')
            end
            %%
            
        end
        
        function obj = rotate(obj)
            %ROTATE Rotates FEATURE object 90 degrees anticlockwise.
            %   Function has no input
            obj.Z = rot90(obj.Z);
        end
        
        function plot(obj)
            %PLOT Displays FEATURE object.
            %   Function has no input
            mesh(obj.X,obj.Y,obj.Z)
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end
end

