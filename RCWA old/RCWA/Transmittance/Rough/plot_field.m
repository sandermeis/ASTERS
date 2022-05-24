function [Er,Eth,Eph] = plot_field(s,beta,labda_x,labda_y,Nx,Ny,P,Q)
X=linspace(0,labda_x,Nx);
Y=linspace(0,labda_y,Ny);

nxc = ceil(Nx/2);
nx1 = nxc - floor(P/2);
nx2 = nxc + floor(P/2);
nyc = ceil(Ny/2);
ny1 = nyc - floor(Q/2);
ny2 = nyc + floor(Q/2);
for i=1:3
sf = zeros(Nx,Ny);
sf(nx1:nx2,ny1:ny2) = s{i};

az = ifft2(ifftshift(sf));
az = az/max(abs(az(:)));
phase = exp(-1i*(beta(1)*X + beta(2)*Y));

E{i} = phase.*az;
end
[Er,Eth,Eph]=plot_spherical(E);
scatter3(Eth(:),Eph(:),Er(:))
% %figure
% % [x,y]=meshgrid(X,Y);
% % quiver3(x,y,0*y,E{1}*conj(E{1}),E{2}*conj(E{2}),E{3}*conj(E{3}))
% %for i=1:3
% figure
% imagesc(X,Y,real(E{1}*conj(E{1})+E{2}*conj(E{2})+E{3}*conj(E{3})))
% ax = gca;
% ax.YDir = 'normal';
% %end
end

