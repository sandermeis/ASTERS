clear all
close all
clc

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS"];

n = ["GaAs"];
L = [300];	% Length vector (nm)
len_n=length(L);

for i=1:length(materials)
    FileName=strcat(["../ri/refractive indices.xlsx - "],materials{i},".csv");
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

eps_r = [eps_rk];
mu_r	= ones(len_n);           % Permeability vector

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
Vr{i}=V;

X=exp(-1i*k_z(i)*k_0*L(i));

S{1,1}=inv(A-X*B*inv(A)*X*B)*(X*B*inv(A)*X*A-B);
S{2,2}=S{1,1};

S{1,2}=inv(A-X*B*inv(A)*X*B)*X*(A-B*inv(A)*B);
S{2,1}=S{1,2};

S_global = RH_star(S_global,S);

Sg{i}=S_global;

end

W_ref=W_0;
W_trn=W_0;

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


c_mn=(Sg{1}{1,2})\(c_ref-Sg{1}{1,1}*c_inc);
c_pl=Sg{1}{2,1}*c_inc+Sg{1}{2,2}*c_mn;
c{1}=0.5*A*c_pl+0.5*B*c_mn;
c{2}=0.5*B*c_pl+0.5*A*c_mn;

s=W*c{1}+W*c{2};
u=-V*c{1}+V*c{2};
lab=1i*k_z*eye(2);
s2=W*expm(-lab*k_0*L)*c{1}+W*expm(lab*k_0*L)*c{2};
u2=-V*expm(-lab*k_0*L)*c{1}+V*expm(lab*k_0*L)*c{2};

J(iter)=real(s(1).*conj(1i*u(2))-s(2).*conj(1i*u(1)))-real(s2(1).*conj(1i*u2(2))-s2(2).*conj(1i*u2(1)));
     

R(iter)=dot(E_ref,E_ref)/dot(E_inc,E_inc);
T(iter)=dot(E_trn,E_trn)/dot(E_inc,E_inc)*real(mu_r(1)/mu_r(end)*k_z(end)/k_z(1));
end

Abs(:,1)=J;
Abs(:,2)=R;
Abs(:,3)=T;
pll=[300:1000';300:1000';300:1000']';
area(pll,Abs)
legend("GaAs","R","T")