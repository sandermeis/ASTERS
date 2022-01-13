clear all
close all
tic
%for roughness=[0.25,0.5,1,2]


%AngleGrid=linspace(-pi/2,pi/2,100);

roughlength=100;
roughdim=10;
roughness=1;

% % No rough layer
% n = ["InGaP","Al03GaAs","InGaP","GaP","GaP","Ag"];
% shape = ["Uniform","Uniform","Uniform","Uniform","GratingXY","Uniform"];
% L = [0,125,50,50,500,50,350,0];

% Rough layer
n = ["InGaP","Al03GaAs","InGaP","GaP","Ag","Ag"];
shape = ["Uniform","Uniform","Uniform","Uniform","Rough","Uniform"];
L = [0,125,50,50,500,roughlength/roughdim*ones(1,roughdim),100,0];

lab1 = 350;
lab2 = 850;
dlab = 2;
lam0_r = lab1:dlab:lab2;

eps_lab = get_lab(lab1,lab2,dlab);

wb = waitbar(0,'Please wait...');

theta	= 0;%8/360*2*pi;
phi     = 0;
pte     = 0.5;
ptm     = 0.5;
rscale = 4;
Xresolution = 10;
Yresolution = 10;
labda_x = 1280*rscale;%120;
labda_y = 1280*rscale;%120;
lenx    = labda_x/Xresolution;
leny    = labda_y/Yresolution;
P       = 9;%5;
Q       = 9;%5;
num_H   = P*Q;
M       = -(P-1)/2:(P-1)/2;
N       = -(Q-1)/2:(Q-1)/2;

% generate rough surface once (not new random every wavelength)

for iter=1:length(lam0_r)%(lam0_r-lab1+1)
    
    k_0 = 2*pi./lam0_r(iter);
    
    [eps_r,mu_r] = build_perm(n,shape,lenx,leny,P,Q,eps_lab,iter,roughdim,roughness);
    [Kx,Ky,Kz] = calc_K(k_0,theta,phi,eps_r,mu_r,M,N,labda_x,labda_y);
    
    s_inc=get_sinc(k_0,eps_r{1},mu_r{1},phi,theta,pte,ptm,num_H);
    
    [r,t,J_im] = calc_layers_transmittance(num_H,mu_r,eps_r,Kx,Ky,Kz,k_0,L,s_inc);
    J(iter,:)=[J_im{:}];
    
    r{3} = -inv(Kz{1})*(Kx*r{1}+Ky*r{2});
    t{3} = -inv(Kz{end})*(Kx*t{1}+Ky*t{2});
    kz0=cos(theta);


    [Rtot(iter),R] = get_power(r,mu_r{1},mu_r{1},Kz{1},kz0);
    [Ttot(iter),T] = get_power(t,mu_r{end},mu_r{1},Kz{end},kz0);
    haze(iter) = (sum(R)-R((num_H+1)/2))/sum(R);

    waitbar(iter/length(lam0_r),wb)
end

toc
delete(wb)

h(1) = plot_results(lam0_r,n,J,Rtot,Ttot);
h(2) = figure;
plot(lam0_r,Rtot)
hold on
plot(lam0_r,haze)
filename="scale4harmonics9roughness"+roughness+".fig";
%savefig(h,filename)
close(h)
%end