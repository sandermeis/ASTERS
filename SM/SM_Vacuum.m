clear all
close all
clc

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS"];

n=["GaAs"];
len_n=length(n);

for i=1:length(materials)
    FileName=strcat(["ri/refractive indices.xlsx - "],materials{i},".csv");
  x{i} = csvread(FileName(1,1));
end

lam0=300:1000;

R=zeros(size(lam0));
T=zeros(size(lam0));

for iter=1:length(lam0)

%% Constants

theta	= 0;                % Angle with respect to normal 0<theta<pi/2
phi     = 0;                % Angle with respect to -x 0<theta<2pi
pte     = 1;                % TE component
ptm     = 0;                % TM component
lam     = lam0(iter);                % Vacuum wavelength (nm)

k_0=2*pi/lam;

% Permittivity vector

for i=1:len_n
  eps_rk(i)=getdata(x,materials,n(i),lam);
end

eps_r = [1,eps_rk,1];
mu_r	= ones(1,len_n+2);           % Permeability vector
L	= [0,300,0];	% Length vector (nm)

%% Simulation

% Starting parameters
k_x = sqrt(eps_r(1)*mu_r(1))*sin(theta)*cos(phi);
k_y = sqrt(eps_r(1)*mu_r(1))*sin(theta)*sin(phi);

S_global = {zeros(2),eye(2);eye(2),zeros(2)};

% Free space
k_z_0 = calc_kz(1,1,k_x,k_y);
Q_0=calc_Q(1,1,k_x,k_y);
Omega_0=1i*k_z_0*eye(2);
V_0=Q_0/Omega_0;
W_0=eye(2);

% Loop through layers
for i=1:length(L)

k_z(i) = calc_kz(eps_r(i),mu_r(i),k_x,k_y);
Q = calc_Q(eps_r(i),mu_r(i),k_x,k_y);
Omega=1i*k_z(i)*eye(2);
W=eye(2);
V=Q*W*inv(Omega);

A=inv(W)*W_0+inv(V)*V_0;
B=inv(W)*W_0-inv(V)*V_0;

Wr{i}=W;
Ar{i}=A;
Br{i}=B;

lam=1i*k_z(i)*k_0*L(i);

X=exp(lam);

S{1,1}=inv(A-X*B*inv(A)*X*B)*(X*B*inv(A)*X*A-B);
S{2,2}=S{1,1};

S{1,2}=inv(A-X*B*inv(A)*X*B)*X*(A-B*inv(A)*B);
S{2,1}=S{1,2};

S_global = RH_star(S_global,S);

end

% %% ref
% 
% k_z_ref = calc_kz(eps_r1,mu_r1,k_x,k_y);
% Q_ref=calc_Q(eps_r1,mu_r1,k_x,k_y);
% 
% Omega_ref=1i*k_z_ref*eye(2);
% 
% V_ref=Q_ref*inv(Omega_ref);
% W_ref=eye(2);
% 
% A_ref=inv(W_0)*W_ref+inv(V_0)*V_ref;
% B_ref=inv(W_0)*W_ref-inv(V_0)*V_ref;
% 
% S_ref{1,1}=-A_ref\B_ref;
% S_ref{1,2}=2*inv(A_ref);
% S_ref{2,1}=0.5*(A_ref-B_ref/A_ref*B_ref);
% S_ref{2,2}=B_ref/A_ref;
% 
% %% trn
% 
% k_z_trn = calc_kz(eps_r2,mu_r2,k_x,k_y);
% Q_trn=calc_Q(eps_r2,mu_r2,k_x,k_y);
% 
% Omega_trn=1i*k_z_trn*eye(2);
% 
% V_trn=Q_trn*inv(Omega_trn);
% 
% W_trn=eye(2);
% 
% A_trn=inv(W_0)*W_trn+inv(V_0)*V_trn;
% B_trn=inv(W_0)*W_trn-inv(V_0)*V_trn;
% 
% S_trn{1,1}=B_trn*inv(A_trn);
% S_trn{1,2}=0.5*(A_trn-B_trn*inv(A_trn)*B_trn);
% S_trn{2,1}=2*inv(A_trn);
% S_trn{2,2}=-inv(A_trn)*B_trn;
% 
% %%
% S_global = RH_star(S_ref,S_global);
% S_global = RH_star(S_global,S_trn);

W_ref=Wr{1};
W_trn=Wr{end};

% Starting parameters
k_inc = k_0*sqrt(eps_r(1)*mu_r(1))*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]';
normal=[0;0;-1];

TE_direction=cross(k_inc,normal)/norm(cross(k_inc,normal));
nancheck=all(isnan(TE_direction));
TE_direction(isnan(TE_direction))=0;
a_TE=TE_direction+nancheck*[0,1,0]';
a_TM=cross(a_TE,k_inc)/norm(cross(a_TE,k_inc));

p=pte*a_TE+ptm*a_TM;

p=p/norm(p);

c_inc=inv(W_ref)*p(1:2);
E_inc=p;

c_ref=S_global{1,1}*c_inc;
c_trn=S_global{2,1}*c_inc;

E_ref=W_ref*c_ref;
E_trn=W_trn*c_trn;

E_ref(3)=-(k_x*E_ref(1)+k_y*E_ref(2))/k_z(1);
E_trn(3)=-(k_x*E_trn(1)+k_y*E_trn(2))/k_z(end);
E_inc(3)=-(k_x*E_inc(1)+k_y*E_inc(2))/k_z(1);

R(iter)=dot(E_ref,E_ref)/dot(E_inc,E_inc);

T(iter)=dot(E_trn,E_trn)/dot(E_inc,E_inc)*real(mu_r(1)/mu_r(end)*k_z(end)/k_z(1));

% % Loop through layers
% 
% c_plus{1}=c_inc;
% c_min{1}=c_ref;
% for i=1:length(L)
% c_plus{i+1}=0.5*(Ar{i}*c_plus{i}*exp(-1i*k_z(i)*k_0*L(i))+Br{i}*c_min{i}*exp(1i*k_z(i)*k_0*L(i)));
% c_min{i+1}=0.5*(Br{i}*c_plus{i}*exp(-1i*k_z(i)*k_0*L(i))+Ar{i}*c_min{i}*exp(1i*k_z(i)*k_0*L(i)));
% 
% end
% 
% E_g=c_plus{2}+c_min{2};
% E_g(3)=-(k_x*E_g(1)+k_y*E_g(2))/calc_kz(eps_r(2),mu_r(2),k_x,k_y);
% 
% G(iter)=dot(E_g,E_g)/dot(E_inc,E_inc);

end

Abs(:,1)=1-R-T;
Abs(:,2)=R;
Abs(:,3)=T;
pll=[300:1000';300:1000';300:1000']';
area(pll,Abs)
legend("GaAs","R","T")