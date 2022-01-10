clear all
close all
tic

AngleGrid=linspace(-pi/2,pi/2,100);

RoughDim=20;

% % No rough layer
% n = ["InGaP","Al03GaAs","InGaP","GaP","Ag"];
% shape = ["Uniform","Uniform","Uniform","Uniform","Uniform"];
% L = [0,125,50,50,500,350,0];

% Rough layer
n = ["InGaP","Al03GaAs","InGaP","GaP","Ag","Ag"];
shape = ["Uniform","Uniform","Uniform","Uniform","Rough","Uniform"];
L = [0,125,50,50,500,5*ones(1,RoughDim),100,0];

lab1 = 350;
lab2 = 850;
lam0_r = lab1:lab2;

eps_lab = get_lab(lab1,lab2);

wb = waitbar(0,'Please wait...');

theta	= 0;
phi     = 0;
pte     = 0.5;
ptm     = 0.5;
labda_x = 10000;
labda_y = 10000;
lenx    = 512;
leny    = 512;
P       = 5;
Q       = 5;
num_H   = P*Q;
M       = -(P-1)/2:(P-1)/2;
N       = -(Q-1)/2:(Q-1)/2;

% generate rough surface once (not new random every wavelength)

for iter=(lam0_r-lab1+1)
    
    k_0 = 2*pi./lam0_r(iter);
    
    [eps_r,mu_r] = build_perm(n,shape,lenx,leny,P,Q,eps_lab,iter);
    [Kx,Ky,Kz] = calc_K(k_0,theta,phi,eps_r,mu_r,M,N,labda_x,labda_y);
    
    s_inc=get_sinc(k_0,eps_r{1},mu_r{1},phi,theta,pte,ptm,num_H);
    
    [r,t,J_im] = calc_layers_transmittance(num_H,mu_r,eps_r,Kx,Ky,Kz,k_0,L,s_inc);
    J(iter,:)=[J_im{:}];
    
    r{3} = -inv(Kz{1})*(Kx*r{1}+Ky*r{2});
    t{3} = -inv(Kz{end})*(Kx*t{1}+Ky*t{2});
    kz0=cos(theta);


    [Rtot(iter),R] = get_power(r,mu_r{1},mu_r{1},Kz{1},kz0);
    [Ttot(iter),T] = get_power(t,mu_r{end},mu_r{1},Kz{end},kz0);
    haze(iter) = sum(R)/R((num_H+1)/2)-1;

    waitbar(iter/length(lam0_r),wb)
end

toc
delete(wb)

plot_results(lam0_r,n,J,Rtot,Ttot)
figure
plot(lam0_r,Rtot)
hold on
plot(lam0_r,haze)