function Zout = epsmatrix(t,eps,xl,yl)

Zin = eps*ones(xl,yl);

if t==1
% Grating
pgon = polyshape([0 0 0.5*yl 0.5*yl], [0 xl xl 0]);
elseif t==2
% Triangle
pgon = polyshape([0.1*yl 0.5*yl 0.9*yl], [0.1*xl 0.9*xl 0.1*xl]);
elseif t==3
%Circle
n=20;
theta = (0:n-1)*(2*pi/n);
r = 0.35*xl;
xc = 0.5*yl;
yc = 0.5*xl;
x = xc + r*cos(theta);
y = yc + r*sin(theta);
pgon = polyshape(x,y);
end
[colGrid,rowGrid] = meshgrid(1:size(Zin,2),1:size(Zin,1));
idx = isinterior(pgon,[colGrid(:),rowGrid(:)]);
idx = reshape(idx,size(Zin));

Zout = Zin;
Zout(idx) = 1;
end

