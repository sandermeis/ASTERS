clear all
close all
tic
%for roughness=2:2:10


%AngleGrid=linspace(-pi/2,pi/2,100);

roughlength=200;
roughdim=20;
roughness=10;

% % No rough layer
% n = ["InGaP","Al03GaAs","InGaP","GaP","GaP","Ag"];
% shape = ["Uniform","Uniform","Uniform","Uniform","GratingXY","Uniform"];
% L = [0,125,50,50,500,50,350,0];

% Rough layer
n = ["InGaP","Al03GaAs","InGaP","GaP","Ag","Ag"];
shape = ["Uniform","Uniform","Uniform","Uniform","Rough","Uniform"];
L = [0,125,50,50,400,roughlength/roughdim*ones(1,roughdim),100,0];

lab1 = 600;
lab2 = 700;
dlab = 10;
lam0_r = lab1:dlab:lab2;

eps_lab = get_lab(lab1,lab2);

wb = waitbar(0,'Please wait...');

theta	= 8/360*2*pi;
phi     = 0;
pte     = 0.5;
ptm     = 0.5;
Xresolution = 10;
Yresolution = 10;
labda_x = 5120;
labda_y = 5120;
lenx    = labda_x/Xresolution;
leny    = labda_y/Yresolution;
P       = 7;
Q       = 7;
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

plot_results(lam0_r,n,J,Rtot,Ttot)
figure
plot(lam0_r,Rtot)
hold on
plot(lam0_r,haze)
filename="Roughness"+roughness+".fig";
%savefig(filename)
%end