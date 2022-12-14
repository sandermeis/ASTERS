function RCWA_plot(param, Sz, layer, sim_num, titlestring, whichdisp)
arguments
param
Sz
layer
sim_num
titlestring = "Sim " + string(sim_num)
whichdisp = 3
end

%n = fixLayerString(layer, param);
%n = {layer.material};
switch whichdisp
    case 1
        plot_results(param, Sz, layer, sim_num, titlestring);
    case 2
        plothaze(param, Sz, layer, sim_num, titlestring);
    case 3
        jsc_harmonics(param, Sz, layer, sim_num, titlestring);
    case 4
        plot_results(param, Sz, layer, sim_num, titlestring);
        plothaze(param, Sz, layer, sim_num, titlestring);
        jsc_harmonics(param, Sz, layer, sim_num, titlestring);
end

end


% figure
% [Er,Eth,Eph] = plot_field({reshape(r{1},P,Q),reshape(r{2},P,Q),reshape(r{3},P,Q)},beta,labda_x,labda_y,lenx,leny,P,Q);
% plot_harmonics(Kx,Ky,Kz{1},r)
% h(2) = tile_harmonics(Sz,string({layer.material}));


function n = fixLayerString(layer,param)
n = string({layer.material});
if param.calcAllRough
    for i = 1:numel(layer)
        
        for j = 1:numel(layer(i).L)
            if numel(layer(i).L)>1
                if i==1
                    n2(j,i) = "Rough "+param.ref_medium+"/"+layer(i).material+" ("+string(j)+") "+string(layer(i).L(j))+" nm";
                elseif i==numel(layer)
                    n2(j,i) = "Rough "+layer(i).material+"/"+param.trn_medium+" ("+string(j)+") "+string(layer(i).L(j))+" nm";
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
    
    shape = cellfun(@(x) isa(x,"Surface"),{layer.input});
    sh=(shape~=0);
    sh_shift = circshift(sh, -1, 2);
    sh_shift(:, end) = 0;
    
    if shape(1)~=0
        n(1)="Rough "+param.ref_medium+"/"+n(1);
        sh(1)=0;
        n(sh)="Rough "+n(sh_shift)+"/"+n(sh);
    else
        n(sh)="Rough "+n(sh_shift)+"/"+n(sh);
    end
    
end
end


function plothaze(param, Sz, layer, sim_num, titlestring)

wavelengthArray = param.wavelengthArray;

% why abs???
centralH = squeeze(abs(Sz((end+1)/2,:,:)));
sumH = squeeze(sum(abs(Sz),1));
diffH = sumH - centralH;
haze = arrayfun(@(a, b) a/b * (b>1e-12),diffH,sumH);

n = {layer.material};
n(end+1) = {"R"};
n(end+1) = {"T"};

for j=1:numel(n)
n{j} = [n{j}{:}];
end

figure
colororder(parula(size(haze,1)))
for i=1:size(haze,1)
    hold on
    plot(wavelengthArray,abs(haze(i,:)).','LineWidth',2)
end
xlabel('Wavelength (nm)')
ylabel('Haze')
legend(n,'location','eastoutside')
title(titlestring,"FontSize",16,"FontWeight",'bold')


%%
% h(2) = figure;
% plot(wavelengthArray,Rtot,'Marker', 'o')
% hold on
% plot(wavelengthArray,haze,'Marker', 'o')

% load('RD.mat','RD');
% plot(RD{2,7}(11:141),RD{2,8}(11:141),'Marker', 'o')
% plot(RD{2,7}(11:141),RD{2,9}(11:141),'Marker', 'o')

% xlabel("Wavelength (nm)")
% ylabel("Absorption (a.u.)")
% legend(["R_{RCWA}","Haze_{RCWA}","R_{measured}","Haze_{measured}"])
% xlim([wavelengthArray(1), wavelengthArray(end)])

end


function plot_results(param, Sz_in, layer, sim_num, titlestring)

Sz = squeeze(sum(Sz_in,1)).';

N = size(Sz,2);

x_grid=repmat(param.wavelengthArray',1,N);

n = {layer.material};
n(end+1) = {"R"};
n(end+1) = {"T"};

for j=1:numel(n)
n{j} = [n{j}{:}];
end

%A(:, [1 2]) = A(:, [2 1]);
gaaspos = find(n=="GaAs"|n=="GaAs_1meV"|n=="GaAs_5meV"|n=="GaAs_10meV"|n=="GaAs_15meV"|n=="GaAs_20meV"|n=="GaAs_25meV", 1);
if ~isempty(gaaspos)
    Sz(:, [1 gaaspos]) = Sz(:, [gaaspos 1]);
    n([1 gaaspos]) = n([gaaspos 1]);
end

for i = 1:N
    jsc(i) = Jsc(Sz(:,i)', param.wavelengthArray);
end

h = figure('Color','w','Position', [100 100 1300 600]);

colors = [76,144,186;43,194,194;244,184,17;222,102,62;255,145,43]/255;
%colors = ["#7f58af","#64c5eb","#e84d8a","#feb326"];
%colors = ["#031f4b","#04396c","#035b96","#6497b1","#b3cde0"];
%colors = jet(4);
colororder(colors);
hx = subplot(1,4,[1,2,3]);
% pos = get(hx,'Position');
% pos(1) = 0.055;
% pos(3) = 0.9;
% set(hx, 'Position', pos)
ar = area(x_grid, Sz, 'LineWidth', 2);%;,'EdgeColor','none');

title(titlestring, "FontSize", 18, "FontWeight", 'bold')
xlim([param.wavelengthArray(1), param.wavelengthArray(end)])
ylim([0, 1])
xlabel("Wavelength (nm)", "FontSize", 16, "FontWeight", 'bold')
ylabel("Absorption (a.u.)", "FontSize", 16, "FontWeight", 'bold')
set(hx,'FontSize',14,'LineWidth', 2)

legstring = pad(n) + sprintf(" Jsc: ") + pad(compose('%0.2f',string(jsc)),'left') + sprintf(" mA/cm^2");
hx2 = subplot(1,4,4);
hx2.Visible = 'off';
[~,legend_h,~,~] = legendflex(ar, cellstr(legstring),'ref', hx2, 'anchor', [1 1], 'buffer', [0 0],'FontName','monospaced');

%hold on
%ar2 = area(x_grid,Sz);

hstyle = {'single','single','cross','single','single','cross'};
hdir = {0,45,0,90,135,45};

for i=1:numel(ar)
h_ind = ceil(i/length(colors))-1;
    if h_ind>0 && h_ind<length(colors)
        j = mod(i-1,length(hstyle))+1;
        hatchfill2(ar(i), hstyle{j},'HatchAngle',hdir{j},'HatchDensity',40,'HatchColor','k','HatchLineWidth',2)
        hatchfill2(legend_h(length(ar)+i), hstyle{j},'HatchAngle',hdir{j},'HatchDensity',40,'HatchColor','k','HatchLineWidth',2)
    end
end

end


function jsc_harmonics(param, Sz, layer, sim_num, titlestring)

n = {layer.material};
n(end+1) = {"R"};
n(end+1) = {"T"};

for j=1:numel(n)
n{j} = [n{j}{:}];
end
    
    num_lay = size(Sz, 2);
    h = figure;
    t = tiledlayout(h,'flow');
    title(t,"Jsc per harmonic per material, " + titlestring,"FontSize",16,"FontWeight",'bold')
    cmin_old = 0;
    cmax_old = 0;
    for i = 1:num_lay
        jsc_empty = zeros(param.P,param.Q);
        h(i) = nexttile;
        jsc_empty(param.tr_ind) = Jsc(squeeze(Sz(:,i,:)), param.wavelengthArray); %sum(squeeze(Sz(:,i,:)),2);
        im = imagesc(jsc_empty);
        title(n(i))
        set(h(i),"Color","none",'YDir','normal','TickLength', [0 0])
        im.AlphaData = jsc_empty~=0;
        xticks(1:param.P)
        yticks(1:param.Q)
        xticklabels(h(i),num2cell(-0.5*(param.P-1):0.5*(param.P-1)))
        yticklabels(h(i),num2cell(-0.5*(param.Q-1):0.5*(param.Q-1)))

        h2(i) = axes(t);
        h2(i).Layout.Tile = h(i).Layout.Tile;
        set(h2(i),"Color","none",'YDir','normal')
        h2(i).XLim = [0,param.P];
        h2(i).YLim = [0,param.Q];
        h2(i).YTickLabel = {};
        h2(i).XTickLabel = {};
        xticks(h2(i), 0:param.P)
        yticks(h2(i), 0:param.Q)
        h2(i).XGrid = 'on';
        h2(i).YGrid = 'on';
        h2(i).GridAlpha = 1;
        j_nocentral = jsc_empty;
        j_nocentral(round(param.P/2),round(param.Q/2))=NaN;
        cmin_new = min(j_nocentral, [], 'all');
        cmax_new = max(j_nocentral, [], 'all');
        cmin_old = max(cmin_old, cmin_new);
        cmax_old = max(cmax_old, cmax_new);
    %alphamap();
    %alpha([0;ones(2,1)])
    %alpha('color');
    end
    set(h, 'Colormap', jet, 'CLim', [cmin_old cmax_old])
    cbh = colorbar(h(end));
    cbh.Layout.Tile = 'east';

end



function h = tile_harmonics4d(Sz, wavelengthArray, n, P, Q, tr_ind)

n = [n, "R", "T"];
h = figure;
tiledlayout(h,'flow','TileSpacing','Compact');
cmax=max(Sz(:));
cmin=0;%min(Sz(:));
a = repair_harmonics(Sz,P,Q,tr_ind);

[X,Y,Z] = ndgrid(-((P-1)/2):((P-1)/2),-((Q-1)/2):((Q-1)/2),wavelengthArray);

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

