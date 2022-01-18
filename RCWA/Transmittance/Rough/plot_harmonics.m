function plot_harmonics(Kx,Ky,Kz,intensity,arrow)
sz=sqrt(size(Kx));
kx=reshape(diag(Kx),sz);
ky=reshape(diag(Ky),sz);
kz=reshape(diag(Kz),sz);
I=reshape(intensity,sz);
%Angles
r=sqrt(kx.^2+ky.^2+kz.^2);
theta=acos(kz./r)*360./(2*pi); %degrees

phi=pi*(kx<0)+pi/2*(kx==0);
kxi=kx;
kxi(kx==0)=1;
phi=(phi+atan(ky./kxi))*360./(2*pi);
figure
scatter(theta(I>0),phi(I>0),50,I(I>0),'filled')
colorbar()
xlabel("Theta (Degrees)");
ylabel("Phi (Degrees)");
set(gca,'clim',[0 0.001])
if arrow
figure
Il=-log(I);
Il(isinf(Il)) = 0;
q=quiver3(0*kx,0*ky,0*kz,real(Il.*kx),real(Il.*ky),real(Il.*kz),0);
q.ShowArrowHead = 'off';
end
end

