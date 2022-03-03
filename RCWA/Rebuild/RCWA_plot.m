function h = RCWA_plot(layer,device,lam0_r,Sz)

n = fixLayerString(layer, device);

h(1) = tile_harmonics4d(Sz, lam0_r, n, device.P, device.Q, device.tr_ind);

h(2) = plot_results(squeeze(sum(Sz,1)).', lam0_r, n);

h(3) = plothaze(Sz, lam0_r, n);

end

% figure
% [Er,Eth,Eph] = plot_field({reshape(r{1},P,Q),reshape(r{2},P,Q),reshape(r{3},P,Q)},beta,labda_x,labda_y,lenx,leny,P,Q);
% plot_harmonics(Kx,Ky,Kz{1},r)
% h(2) = tile_harmonics(Sz,string({layer.material}));


function n = fixLayerString(layer,device)
n = string({layer.material});
if device.calcAllRough
    for i = 1:numel(layer)
        
        for j = 1:numel(layer(i).L)
            if numel(layer(i).L)>1
                if i==1
                    n2(j,i) = "Rough "+device.ref_medium+"/"+layer(i).material+" ("+string(j)+") "+string(layer(i).L(j))+" nm";
                elseif i==numel(layer)
                    n2(j,i) = "Rough "+layer(i).material+"/"+device.trn_medium+" ("+string(j)+") "+string(layer(i).L(j))+" nm";
                else
                    n2(j,i) = "Rough "+layer(i-1).material+"/"+layer(i).material+" ("+string(j)+") "+string(layer(i).L(j))+" nm";
                end
            else
                n2(j,i) = layer(i).material+" "+string(layer(i).L(j))+" nm";
            end
        end
    end
    n = rmmissing(n2(:));
else
    
    shape=[layer.shape];
    
    sh=(shape~=0);
    sh_shift = circshift(sh, -1, 2);
    sh_shift(:, end) = 0;
    
    if shape(1)~=0
        n(1)="Rough "+device.ref_medium+"/"+n(1);
        sh(1)=0;
        n(sh)="Rough "+n(sh_shift)+"/"+n(sh);
    else
        n(sh)="Rough "+n(sh_shift)+"/"+n(sh);
    end
    
end
end


function h = plothaze(Sz, lam0_r, n)

% why abs???
centralH = squeeze(abs(Sz((end+1)/2,:,:)));
sumH = squeeze(sum(abs(Sz),1));
diffH = sumH - centralH;
haze = arrayfun(@(a, b) a/b * (b>1e-12),diffH,sumH);

n = [n, "R", "T"];
h = figure;
colororder(hsv(8))
for i=1:size(haze,1)
    hold on
    plot(lam0_r,abs(haze(i,:)).','LineWidth',2)
end
xlabel('Wavelength (nm)')
ylabel('Haze')
legend(n,'location','eastoutside')


%%
% h(2) = figure;
% plot(lam0_r,Rtot,'Marker', 'o')
% hold on
% plot(lam0_r,haze,'Marker', 'o')

% load('RD.mat','RD');
% plot(RD{2,7}(11:141),RD{2,8}(11:141),'Marker', 'o')
% plot(RD{2,7}(11:141),RD{2,9}(11:141),'Marker', 'o')

% xlabel("Wavelength (nm)")
% ylabel("Absorption (a.u.)")
% legend(["R_{RCWA}","Haze_{RCWA}","R_{measured}","Haze_{measured}"])
% xlim([lam0_r(1), lam0_r(end)])

end


function h = plot_results(Sz, lam0_r, n)

h = figure;
lab1 = lam0_r(1);
lab2 = lam0_r(end);

N=size(Sz,2);

x_grid=repmat(lam0_r',1,N);

n(end+1)="R";
n(end+1)="T";

colororder(hsv(8))

area(x_grid,Sz,'EdgeColor','none')
% hold on
% plot(lam0_r,real(eps_lab{1})/max(real(eps_lab{1})))
% plot(lam0_r,imag(eps_lab{1})/max(real(eps_lab{1})))
legend(n,'Location','eastoutside')
%ylim([-0.5,1.5])
xlim([lab1-0.1*(lab2-lab1),lab2+0.1*(lab2-lab1)])
end


function h = tile_harmonics4d(Sz, lam0_r, n, P, Q, tr_ind)

n = [n, "R", "T"];
h = figure;
tiledlayout(h,'flow','TileSpacing','Compact');
cmax=max(Sz(:));
cmin=0;%min(Sz(:));
a = repair_harmonics(Sz,P,Q,tr_ind);

[X,Y,Z] = ndgrid(-((P-1)/2):((P-1)/2),-((Q-1)/2):((Q-1)/2),lam0_r);

for i=1:size(Sz,2)

D=squeeze(a(:,:,i,:));
nexttile
scatter3(Z(:),X(:),Y(:),400,D(:),'filled');
xlabel("Wavelength (nm)")
ylabel("P")
zlabel("Q")
    title(n(i))
    axis ij
    caxis manual
    caxis([cmin cmax]);
     alpha color
     alpha scaled

end
 cb = colorbar;
 cb.Layout.Tile = 'east';

end


function a = repair_harmonics(Sz,P,Q,tr_ind)

num_lay = size(Sz,2);
num_lab = size(Sz,3);
a = zeros(P,Q,num_lay,num_lab);

for i=1:num_lay
    for j=1:num_lab
        c=a(:,:,i,j);
        c(tr_ind)=Sz(:,i,j);
        a(:,:,i,j) = c;
    end
end
end








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


function h = tile_harmonics(Sz,n)

n = [n, "R", "T"];
h = tiledlayout('flow','TileSpacing','Compact');

for i=1:size(Sz,2)
    nexttile
    imagesc(squeeze(Sz(:,i,:)))
    xlabel("Wavelength (nm)")
    ylabel("Harmonics")
    title(n(i))
end

end


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


function thbin = azimuthal_bin(TH,R,th_range)

for i = 1:length(th_range)-1
    thbin(i) = mean(R(TH>=th_range(i) & TH<th_range(i+1)));  
end
end

