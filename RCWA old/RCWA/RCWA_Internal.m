clear all
close all
tic

n = ["GaAs","Ag"];
shape = ["Uniform","Uniform"];
L = [0,300,100,0];

%n=["MgF2","ZnS","Al03GaAs","GaAs","InGaP","Al03GaAs","Ag"];
%shape=["Uniform","Uniform","Uniform","Uniform","Uniform","Uniform","Uniform"];

% MgF2      94
% ZnS       44
% Al03GaAs	25
% GaAs      300
% InGaP     100
% Al03GaAs	100
% Ag        100

lab1 = 300;
lab2 = 900;
lam0_r = lab1:lab2;

eps_lab = get_lab(lab1,lab2);

% for l=1:length(lam0_r)
%     eps_lab{1}(l)=real(eps_lab{1}(l));
% end

wb = waitbar(0,'Please wait...');

theta	= 0;
phi     = 0;
pte     = 0.5;
ptm     = 0.5;
labda_x = 100;
labda_y = 100;
lenx    = 512;
leny    = 512;
P       = 3;
Q       = 3;
num_H   = P*Q;
M       = -(P-1)/2:(P-1)/2;
N       = -(Q-1)/2:(Q-1)/2;

for iter=(lam0_r-lab1+1)
    
    k_0 = 2*pi./lam0_r(iter);
    
    [eps_r,mu_r]=build_perm(n,shape,lenx,leny,P,Q,eps_lab,iter);
    
    %% Simulation
    
    % Starting parameters
    [Kx,Ky,Kz] = calc_K(k_0,theta,phi,eps_r,mu_r,M,N,labda_x,labda_y);
    S_global = {zeros(2*num_H),eye(2*num_H);eye(2*num_H),zeros(2*num_H)};
    
    % Free space
    P_0 = calc_PQ(eye(size(Kx)),eye(size(Kx)),Kx,Ky);
    Q_0 = calc_PQ(eye(size(Kx)),eye(size(Kx)),Kx,Ky);
    [W_0, eigval0]=eig(P_0*Q_0);
    V_0=Q_0*W_0*inv(sqrt(eigval0));
    
    % Loop through layers
    for i=1:length(L)
        [S,W{i},V{i},l{i}] = calc_layer(V_0,W_0,mu_r{i},eps_r{i},Kx,Ky,k_0*L(i));
        S_global = RH_star(S_global,S);
        Sg{i} = S_global;
    end

    [E_inc{1},E_inc{2}] = get_Einc(phi,theta,pte,ptm,num_H);
    E_inc{3} = -inv(Kz{1})*(Kx*E_inc{1}+Ky*E_inc{2});
    
    c_inc=inv(W{1})*[E_inc{1};E_inc{2}];
    c_ref=S_global{1,1}*c_inc;
    c_trn=S_global{2,1}*c_inc;

    [H_inc{1},H_inc{2}]=split_xy(-V{1}*c_inc);
    s_inc=0.5*real(E_inc{1}.*conj(1i*H_inc{2})-E_inc{2}.*conj(1i*H_inc{1}));
    
    
    % Calculate fields by backward propagation
     for i=(length(L)-1):-1:2
        c = get_mode(c_ref,c_inc,Sg{i-1},W_0,V_0,W{i},V{i});
        [s_f,s_b,u_f,u_b] = get_total_field(Kx,Ky,eps_r{i},mu_r{i},W{i},V{i},l{i},c{1},c{2},0);
        [s_f2,s_b2,u_f2,u_b2] = get_total_field(Kx,Ky,eps_r{i},mu_r{i},W{i},V{i},l{i},c{1},c{2},k_0*L(i));
        J{i-1}(iter)=(sum(Sz(sum_elements(s_f,s_b),sum_elements(u_f,u_b)))-sum(Sz(sum_elements(s_f2,s_b2),sum_elements(u_f2,u_b2))))/sum(s_inc);
     end
    
    E_ref=get_field(Kx,Ky,Kz{1},W{1},0,c_ref,0*c_ref,0);
    E_trn=get_field(Kx,Ky,Kz{end},W{end},0,c_trn,0*c_trn,0);

    Rtot(iter) = get_power(E_ref,mu_r{1},mu_r{1},Kz{1},cos(theta));
    Ttot(iter) = get_power(E_trn,mu_r{end},mu_r{1},Kz{end},cos(theta));

    waitbar(iter/length(lam0_r),wb)
end
toc

delete(wb)
%J{1}=0*J{1};

plot_results(lam0_r,n,J,Rtot,Ttot)