clear all
close all
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
%%% 7 = Real surface, only 512x512, 10umx10um

n           = ["InGaP","Al03GaAs","InGaP","GaP","Ag","Ag"];
shape       = [0,0,0,0,7,0];
roughdim	= [0,0,0,0,25,0];
L           = [125,50,50,125,425,100];

roughness = 1;

lab1 = 650;
lab2 = 655;
dlab = 5;
lam0_r = lab1:dlab:lab2;

eps_lab = get_lab(lab1,lab2,dlab);

wb = waitbar(0,'Please wait...');

theta	= 0;%82/360*2*pi;%pi/4;
phi     = 0;
pte     = 0.5;
ptm     = 0.5;
% rscale = 4;
% Xresolution = 10;
% Yresolution = 10;
labda_x = 10000;%1280*rscale;%120;
labda_y = 10000;%1280*rscale;%120;
lenx    = 512;%labda_x/Xresolution;
leny    = 512;%labda_y/Yresolution;
P       = 9;
Q       = 9;
num_H   = P*Q;

for iter=1:length(lam0_r)
    
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

    waitbar(iter/length(lam0_r),wb)
end

toc
delete(wb)

% figure
% [Er,Eth,Eph]=plot_field({reshape(r{1},P,Q),reshape(r{2},P,Q),reshape(r{3},P,Q)},beta,labda_x,labda_y,lenx,leny,P,Q);

plot_harmonics(Kx,Ky,Kz{1},R,true)

h(1) = plot_results(lam0_r,n,shape,J,Rtot,Ttot);
h(2) = figure;
plot(lam0_r,Rtot)
hold on
plot(lam0_r,haze)
%filename="scale4harmonics9roughness"+roughness+".fig";
%savefig(h,filename)
%close(h)
