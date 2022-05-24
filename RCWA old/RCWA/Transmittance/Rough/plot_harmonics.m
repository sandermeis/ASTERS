function plot_harmonics(Kx,Ky,Kz,r)

sz=sqrt(size(Kx));
kx=reshape(diag(Kx),sz);
ky=reshape(diag(Ky),sz);
kz=reshape(diag(Kz),sz);

I=r{1}.*conj(r{1})+r{2}.*conj(r{2})+r{3}.*conj(r{3});
Sx=reshape(real(Kx*I),sz);
Sy=reshape(real(Ky*I),sz);
Sz=reshape(real(Kz*I),sz);

% integral(S(r,th,ph)*r*sin(theta)dth dr)

%Angles
% figure
% scatter(Sx(:),Sy(:),Sz(:),50,'filled')
% colorbar
% figure
% q=quiver3(0*Sx(:),0*Sy(:),0*Sz(:),Sx(:),Sy(:),Sz(:),0);
% q.ShowArrowHead = 'off';
% figure
% q=quiver3(0*Kx(:),0*Ky(:),0*Kz(:),Kx(:),Ky(:),Kz(:),0);
% q.ShowArrowHead = 'off';
[kph,kth,kr] = cart2sph(kx,ky,kz);
kth=pi/2-kth; %from z to -z
figure
scatter(kth(:),kph(:),50,kr(:),'filled')
title("K")
set(gca,'XTick',0:pi/16:pi/2) 
set(gca,'XTickLabel',{'0','pi/16','pi/8','3pi/16','pi/4','5pi/16','3pi/8','7pi/16','pi/2'})
set(gca,'YTick',-pi:pi/2:pi) 
set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
colorbar
[Sph,Sth,Sr] = cart2sph(Sx,Sy,Sz);
Sth=pi/2-Sth; %from z to -z
figure
scatter(Sth(:),Sph(:),50,Sr(:),'filled')
title("S")
set(gca,'XTick',0:pi/16:pi/2) 
set(gca,'XTickLabel',{'0','pi/16','pi/8','3pi/16','pi/4','5pi/16','3pi/8','7pi/16','pi/2'})
set(gca,'YTick',-pi:pi/2:pi) 
set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
xlabel("Theta")
%xlim([0, pi])
ylabel("Phi")
ylim([-pi, pi])
colorbar()
% figure
% binscatter(Sth(:),Sph(:))
figure
thrange=0:pi/128:pi/2;
thbin=azimuthal_bin(Sth,Sr,thrange);
thbin(isnan(thbin))=0;
plot(thrange(2:end),thbin)

% F = scatteredInterpolant(Sth(:),Sph(:),Sr(:),'nearest','nearest');
% 
% xgrid=0:pi/400:pi/2;
% ygrid=-pi:pi/100:pi;
% [Xg,Yg]=meshgrid(xgrid,ygrid);
% figure
% ZZ=F(Xg,Yg);
% h2=surf(Xg,Yg,ZZ);
% h2.EdgeColor='none';
% zlim([0,0.0005])

% r=sqrt(kx.^2+ky.^2+kz.^2);
% theta=acos(kz./r)*360./(2*pi); %degrees
% 
% phi=pi*(kx<0)+pi/2*(kx==0);
% kxi=kx;
% kxi(kx==0)=1;
% phi=(phi+atan(ky./kxi))*360./(2*pi);
% figure
% scatter(theta(I>0),phi(I>0),50,I(I>0),'filled')
% colorbar()
% xlabel("Theta (Degrees)");
% ylabel("Phi (Degrees)");
% set(gca,'clim',[0 0.001])
% % if arrow
% % figure
% % Il=-log(I);
% % Il(isinf(Il)) = 0;
% % q=quiver3(0*kx,0*ky,0*kz,real(Il.*kx),real(Il.*ky),real(Il.*kz),0);
% % q.ShowArrowHead = 'off';
% % end
end

