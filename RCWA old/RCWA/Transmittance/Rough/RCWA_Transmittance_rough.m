clear all
close all
for gap_length=360:1:360
tic

%%% SHAPE
%%%
%%% 0 = Uniform
%%% 1 = GratingX
%%% 2 = GratingY
%%% 3 = GratingXY
%%% 4 = Triangle
%%% 5 = Circle
%%% 6 = Rough
%%% 7 = Real surface, only 512x512, 10x10 or 30umx30um

% LayerObjects
% Normal
% n=["MgF2", "ZnS", "AlInP", "GaAs", "InGaP", "Al03GaAs", "Ag"];
% shape=[0, 0, 0, 0, 0, 0, 0];
% roughdim=[0, 0, 0, 0, 0, 0, 0];
% L=[94, 44, 20, 300, 100, 100, 100];

% inverse cell
% n=["GaP", "InGaP", "Al03GaAs", "InGaP", "GaAs"];
% shape=[7, 0, 0, 0, 0];
% roughdim=[25, 0, 0, 0, 0];
% L=[400, 50, 50, 125, 3000];

n           = ["GaAs","GaAs","Ag"];
shape       = [0,0,3];
roughdim	= [0,1,0];
L           = [300,100,3000];

roughness = 2;

lab1 = 350;
lab2 = 900;
dlab = 5;
lam0_r = lab1:dlab:lab2;

eps_lab = get_lab(lam0_r);

wb = waitbar(0,'Please wait...');

theta	= 0;%8/360*2*pi;%pi/4;
phi     = 0;
pte     = 0.5;
ptm     = 0.5;
% rscale = 4;
% Xresolution = 10;
% Yresolution = 10;
labda_x = 30000;%1280*rscale;%120;
labda_y = 30000;%1280*rscale;%120;
lenx    = 512;%labda_x/Xresolution;
leny    = 512;%labda_y/Yresolution;
P       = 9;
Q       = 9;
num_H   = P*Q;

for iter=1:length(lam0_r)
    toc_start=toc;
    k_0 = 2*pi./lam0_r(iter);
    
    [eps_r,mu_r] = build_perm(n,shape,lenx,leny,P,Q,eps_lab,iter,roughdim,roughness);
    [Kx,Ky,Kz,beta] = calc_K(k_0,theta,phi,eps_r,mu_r,P,Q,labda_x,labda_y);
    
    s_inc=get_sinc(k_0,eps_r{1},mu_r{1},phi,theta,pte,ptm,num_H);
    
    [r,t,J_im] = calc_layers_transmittance(num_H,mu_r,eps_r,Kx,Ky,Kz,k_0,s_inc,L,shape,roughdim);
    J(iter,:)=[J_im{:}];
    
    r{3} = -inv(Kz{1})*(Kx*r{1}+Ky*r{2});
    t{3} = -inv(Kz{end})*(Kx*t{1}+Ky*t{2});

    [Rtot(iter),R] = get_power(r,mu_r{1},mu_r{1},Kz{1},beta(3)/k_0);
    [Ttot(iter),T] = get_power(t,mu_r{end},mu_r{1},Kz{end},beta(3)/k_0);
    haze(iter) = (sum(R)-R((num_H+1)/2))/sum(R);
    
    toc_arr(iter)=(toc-toc_start);
    timeleft=duration(seconds((length(lam0_r)-iter)*mean(toc_arr)));
    waitbar(iter/length(lam0_r),wb,'Time remaining: ' + string(timeleft,'hh:mm:ss') )
end

toc
delete(wb)

% % figure
% % [Er,Eth,Eph]=plot_field({reshape(r{1},P,Q),reshape(r{2},P,Q),reshape(r{3},P,Q)},beta,labda_x,labda_y,lenx,leny,P,Q);
% %%
% plot_harmonics(Kx,Ky,Kz{1},r)
% %
 h(1) = plot_results(lam0_r,n,shape,J,Rtot,Ttot);
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
%%
% save("haze650ref.mat","haze")
% filename="RCWA_ARough2_"+num_H+"_GaPthickness_"+gap_length+".fig";
% savefig(h,filename)
% clearvars
% close all
end