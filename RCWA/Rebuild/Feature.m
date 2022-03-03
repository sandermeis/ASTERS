classdef Feature
    properties
        size (1,1) {mustBeNumeric} = 625
        resolution = 32 %pts/nm
        height (1,1) {mustBeNumeric} = 16
        shape string = "Pyramid"
        
    end
    properties (SetAccess = private)
        X
        Y
        Z
    end
    
    methods
        function obj = Feature(varargin)
            if nargin==0
                warning('No shape added, producing feature with default parameters')
            elseif nargin>=1
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
                        obj.Z = obj.Z - min(obj.Z(:));
                        obj.resolution  = N(1);
 
                        % rescaling
                        if nargin==2
                            obj.size = varargin{2};
                        elseif nargin==3
                            obj.size = varargin{2};
                            obj.height = varargin{3};
                            obj.Z = obj.height * obj.Z/max(obj.Z(:));
                            warning('Rescaling feature height')
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
            
            switch obj.shape
                case "GratingX"
                    % 1D Grating X
                    obj.Z = obj.height * (obj.X<(obj.size/2));
                case "GratingY"
                    % 1D Grating Y
                    obj.Z = obj.height * (obj.Y<(obj.size/2));
                case "GratingXY"
                    % 2D Grating XY
                    obj.Z = obj.height * ((obj.X<(obj.size/2))&(obj.Y<(obj.size/2)));
                case "Triangle"
                    % 2D Triangle
                    obj.Z = obj.height * ((2*abs(obj.X-(obj.size/2))+obj.Y)<obj.size);
                case "Circle"
                    % 2D Circle
                    obj.Z = obj.height * ((obj.X-obj.size/2).^2+(obj.Y-obj.size/2).^2<=(obj.size/2).^2);
                case "Sphere"
                    % 3D Half sphere
                    tmp = (1-((obj.X-obj.size/2).^2+(obj.Y-obj.size/2).^2)/(obj.size/2)^2);
                    obj.Z = obj.height * (tmp+abs(tmp))/2;
                case "Pyramid"
                    % 3D Pyramid
                    blaat = @(a,b) max([abs(a-obj.size/2),abs(b-obj.size/2)]);
                    obj.Z = obj.height * (1-2*arrayfun(blaat,obj.X,obj.Y)/obj.size);
                case "RidgeX"
                    % 3D Ridge X
                    
                case "RidgeY"
                    % 3D Ridge Y
                case "Custom"
                    
                otherwise
                    warning('No shape added, producing feature with default shape')
            end
            %%
            
        end
        
        function obj = rotate(obj)
            obj.Z = rot90(obj.Z);
        end
        
        function plot(obj)
            mesh(obj.X,obj.Y,obj.Z)
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end
end

