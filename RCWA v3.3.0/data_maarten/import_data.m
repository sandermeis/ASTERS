% im_650=importdata("650C_totaldiffusereflectance_inverted_noARC_QE.txt");
% im_650_ARC=importdata("650C_totaldiffusereflectance_inverted_ARC_centre_QE.txt");
% im_650_AG=importdata("650C_as-grown.xlsx");
% 
% im_700=importdata("700C_totaldiffusereflectance_inverted_noARC_QE.txt");
% im_700_ARC=importdata("700C_totaldiffusereflectance_inverted_ARC_QE.txt");
% im_700_AG=importdata("700C_as-grown.xlsx");
% 
% im_730=importdata("730C_totaldiffusereflectance_inverted_noARC_QE.txt");
% im_730_ARC=importdata("730C_totaldiffusereflectance_inverted_ARC_QE.txt");
% im_730_AG=importdata("730C_as-grown.xlsx");
% 
% im_650_STM=1e9*importdata("650C_aftergrowth_centre_30um.txt");
% im_700_STM=1e9*importdata("700C_aftergrowth_centre_30um.txt");
% im_730_STM=1e9*importdata("730C_aftergrowth_centre_30um.txt");
% 
% im_700_STM_fixed=im_700_STM;
% im_700_STM_fixed(1:118,:)=im_700_STM_fixed(1:118,:)-165;
% im_700_STM_fixed(119,1:232)=im_700_STM_fixed(119,1:232)-165;

%still some funky stuff going on
%plotfig(im_650,im_700,im_730)
%%
im_oxide = importdata("je487-s14-middle-6-OXIDE-SURFACE.txt");
im_algaas = importdata("je484-s14-middle-2-ALGAAS SURFACE.txt");

dat_oxide = im_oxide.data;
dat_algaas = im_algaas.data;

dat_oxide = dat_oxide*1e9;

dat_oxide = [dat_oxide,dat_oxide;dat_oxide,dat_oxide];

dat_algaas = dat_algaas*1e9;
%%
[X,Y]=meshgrid(linspace(0,10000,512),linspace(0,10000,512));

%%
tiledlayout('flow')
nexttile
imagesc(dat_oxide)
colorbar
nexttile
surf(X,Y,dat_algaas,'EdgeColor','none')
axis equal
colorbar
sz = size(dat_oxide);
xg = 1:sz(1);
yg = 1:sz(2);
F = griddedInterpolant({xg,yg},dat_oxide);
xq = (1:2:sz(1))';
yq = (1:2:sz(2))';
dat_oxide_new = F({xq,yq});
nexttile
surf(X,Y,dat_oxide_new,'EdgeColor','none')
axis equal

%writematrix(dat_algaas,'AlGaAs_surface1_10um.csv')
%writematrix(dat_oxide_new,'Oxide_surface1_10um.csv')
%%
% tiledlayout('flow')
% rng(1337)
% N = [512, 512];
% H=1;
% [X,Y] = ndgrid(1:N(1),1:N(2));
% i = min(X-1,N(1)-(X-1));
% j = min(Y-1,N(2)-(Y-1));
% 
% H=0.25;
% Z1 = real(ifft2(exp(-0.5*(i.^2+j.^2)/H^2).*fft2(randn(N))));
% H =2;
% Z2 = real(ifft2(exp(-0.5*(i.^2+j.^2)/H^2).*fft2(randn(N))));
% H =10;
% Z3 = real(ifft2(exp(-0.5*(i.^2+j.^2)/H^2).*fft2(randn(N))));
% nexttile
% [q , C, PSD] = psd_2D(Z1 , 1);
% loglog(q,C)
% nexttile
% [q , C, PSD] = psd_2D(Z2 , 1);
% loglog(q,C)
% nexttile
% [q , C, PSD] = psd_2D(Z3 , 1);
% loglog(q,C)
% 

%%
% RD{1,1}=im_650.textdata(2,3);
% RD{1,2}=im_650.textdata(2,4);
% RD{1,3}=im_650.textdata(2,6);
% RD{1,4}=im_700.textdata(2,1);
% RD{1,5}=im_700.textdata(2,2);
% RD{1,6}=im_700.textdata(2,4);
% RD{1,7}=im_730.textdata(2,3);
% RD{1,8}=im_730.textdata(2,4);
% RD{1,9}=im_730.textdata(2,6);
% 
% RD{2,1}=im_650.data(1:end,3);
% RD{2,2}=im_650.data(1:end,4);
% RD{2,3}=im_650.data(1:end,6);
% RD{2,4}=im_700.data(1:end,1);
% RD{2,5}=im_700.data(1:end,2);
% RD{2,6}=im_700.data(1:end,4);
% RD{2,7}=im_730.data(1:end,3);
% RD{2,8}=im_730.data(1:end,4);
% RD{2,9}=im_730.data(1:end,6);
% 
% RD{3,1}=im_650_STM;
% RD{3,4}=im_700_STM_fixed;
% RD{3,7}=im_730_STM;
% 
% save('RD.mat','RD')
% fsc=importdata("TMM_res.xlsx");
% area(fsc.data.FullSolarCell(1:701,1),fsc.data.FullSolarCell(1:701,2:10));
% legend(fsc.textdata.FullSolarCell(1,2:10))

% fscf=importdata("TMM_res_fix.xlsx");
% figure;
% area(fscf.data(1:701,1),fscf.data(1:701,2:10));
% legend(fscf.textdata(2,2:10))

% load("A2.mat","A2","labda_range","layer_names");
% figure
% area(labda_range,A2);
% title("TMM Absorption")
% xlabel("Wavelength (nm)")
% ylabel("Absorption (a.u.)")
% legend([layer_names{2:8},"R","T"])
% 
% figure
% area(labda_range,[J,Rtot',Ttot']);
% title("RCWA Absorption")
% xlabel("Wavelength (nm)")
% ylabel("Absorption (a.u.)")
% legend([n,"R","T"])
% 
% figure
% imagesc(([J,Rtot',Ttot']-A2));
% c = colorbar;  
% % c.Ruler.TickLabelFormat='%.0e';
% % c.Ruler.Exponent = 0;
% title("Difference in TMM and RCWA Absorption")
% kz=1:10:length(labda_range);
% yticks(kz)
% yticklabels(labda_range(kz))
% xticklabels([n,"R","T"])
% ylabel("Wavelength (nm)")

%% FEATURE 1 DATA PROCESSSING
a=im_700_STM_fixed;
%a(im_650_STM>100)=NaN;
t=tiledlayout('flow');
nexttile
mesh(linspace(0,10000,512),linspace(0,10000,512),a)
%grid=linspace(10,5120,512);
% mesh(grid(200:275),grid(256:382),a(256:382,200:275))
% axis equal
%writematrix(a(256:382,200:275),'feature1.csv')

nexttile
b=a(256:382,200:275);
mesh(b)
nexttile
b=b-min(b(:));
b(b<200)=NaN;
rowgrid = linspace(0,10000/512*size(b,1),size(b,1));
colgrid = linspace(0,10000/512*size(b,2),size(b,2));
[X,Y] = meshgrid(rowgrid,colgrid);
X=X';
Y=Y';
mesh(X,Y,b)
nexttile
%c=abs(gradient(abs(gradient(b))));
%contourf(b,30,'ShowText','on')

k = boundary(X(:),Y(:));
[kx,ky]=ind2sub(size(b),k);
% X(k) Y(k)
u=find(~isnan(b));
[ux,uy]=ind2sub(size(b),u);
% [X,Y] = meshgrid(-3:3)
% [Xq,Yq] = meshgrid(-3:0.25:3)
Xk=[X(k);X(u)];
Yk=[Y(k);Y(u)];
bk=[zeros(size(b(k)));b(u)];
V = scatteredInterpolant(Xk,Yk,bk);
Z=V(X,Y);
Z(Z<0)=0;
scatter3(Xk,Yk,bk)
nexttile
mesh(X,Y,Z)
% pad zeros
Z2=zeros(128,128);
Z2(1:127,27:102)=Z;
%writematrix(Z2,'Feature1_2500nm.csv')
%%
function plotfig(varargin)
switch nargin
    case 1
        im=varargin{1};
        if ~isstruct(varargin{1})
            
            figure
            imagesc(im*1e9)
            colorbar()
        else
            figure
            leg=[];
            for i=1:(length(im.colheaders))
                if im.colheaders{i}=="Wavelength [nm]"
                    plot(im.data(:,i),im.data(:,i+1));
                    if isempty(leg)
                        leg=[string(im.colheaders{i+1})];
                    else
                        leg=[leg,string(im.colheaders{i+1})];
                    end
                    hold on
                end
            end
            legend(leg)
        end
    otherwise
        figure
        leg=[];
        hold on
        for j=1:nargin
            im=varargin{j};
            kz=3-2*(j==2);
            for i=kz:kz%(length(im.colheaders))
                if im.colheaders{i}=="Wavelength [nm]"
                    plot(im.data(:,i),im.data(:,i+3)./im.data(:,i+1),'Marker', 'o');
                    if isempty(leg)
                        leg=[string(im.colheaders{i+1})];
                    else
                        leg=[leg,string(im.colheaders{i+1})];
                    end
                    
                end
            end
            
        end
        legend(leg)
end
end
