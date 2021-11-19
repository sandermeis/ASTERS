function Zout = epsmatrix(xl,yl)
% 100-by-100 sample matrix
Zin = ones(xl,yl);
% Example of polygon
pgon = polyshape([0.1*yl 0.5*yl 0.9*yl], [0.1*xl 0.9*xl 0.1*xl]);
% Create logical index which indicates inside/outsize the polygon
[colGrid,rowGrid] = meshgrid(1:size(Zin,2),1:size(Zin,1));
idx = isinterior(pgon,[colGrid(:),rowGrid(:)]);
idx = reshape(idx,size(Zin));
% Set elements inside the polygon to 0
Zout = Zin;
Zout(idx) = 0;
end

