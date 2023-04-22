% rm=readmatrix("data_daan/730C_aftergrowth_centre_30um.txt");
% rm=rm*1e9;
% 
% maxrm=max(rm,[],'all');
% 
% [dnew, Lnew, Lrecalc]=discretize_surface(rm, 64, 0, false, false, 0, 0);
% dnew(dnew==2)=0;

% ADD PBC
function z2 = gauss2d(N, mux_in, muy_in, sx, sy, h)
%[Nx,Ny,Nz]=size(dnew);

% sx=Nx/8;
% sy=Nx/8;
rho=0;
[x,y]=ndgrid(1:N,1:N);
% mux=Nx/2;
% muy=Nx/2;
% h=maxrm/5;
%%
mux=round(N/2);
muy=round(N/2);
z0=floor(h*exp(-1/(2*(1-rho^2))*(((x-mux)/sx).^2+((y-muy)/sy).^2-2*rho*(x-mux)/sx*(y-muy)/sy)));
z1=circshift(z0,mux_in-mux,1);
z2=circshift(z1,muy_in-muy,2);
end
% tiledlayout('flow')
% nexttile
% surf(x,y,rm,'EdgeColor','none')
% nexttile
% surf(x,y,z,'EdgeColor','none')
% nexttile
% surf(x,y,rm-z,'EdgeColor','none');