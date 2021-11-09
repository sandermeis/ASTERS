clear all
close all
clc

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS"];

n=["GaAs"];
len_n=length(n);

for i0=1:length(materials)
    FileName=strcat(["refractive indices.xlsx - "],materials{i0},".csv");
  x{i0} = csvread(FileName(1,1));
end

lam0=300:1000;

R=zeros(size(lam0));
T=zeros(size(lam0));

for iter=1:length(lam0)
% Abs,          ref,             trn
%0.5225478705	0.4774521289	5.57E-10
%% Constants

eps_r1     = 1+0i;                % Permittivity reflection medium
eps_r2     = 1+0i;                % Permittivity transmission medium
mu_r1     = 1+0i;                % Permeability reflection medium
mu_r2     = 1+0i;                % Permeability transmission medium
theta	= 0;                % Angle with respect to normal 0<theta<pi/2
phi     = pi/4;                % Angle with respect to -x 0<theta<2pi
pte     = 1;                % TE component
ptm     = 0;                % TM component
lam     = lam0(iter);                % Vacuum wavelength (nm)

k_0=2*pi/lam;

for i1=1:len_n
  eps_r(i1)=getdata(x,materials,n(i1),lam);
end
%eps_r	= [13.711+18.298*1i];      % Permittivity vector
mu_r	= ones(1,len_n);           % Permeability vector
L	= [300];	% Length vector (nm)

%% Simulation

% Starting parameters
k_x = sqrt(eps_r1*mu_r1)*sin(theta)*cos(phi);
k_y = sqrt(eps_r1*mu_r1)*sin(theta)*sin(phi);

S_global = {zeros(2),eye(2);eye(2),zeros(2)};

% add boundaries
eps_rk=[1,eps_r,1];
mu_rk=[1,mu_r,1];

% Loop through layers
for i2=1:length(eps_rk)

k_z(i2) = calc_kz(eps_rk(i2),mu_rk(i2),k_x,k_y);
Q=calc_Q(eps_rk(i2),mu_rk(i2),k_x,k_y);
V{i2}=Q*inv(1i*k_z(i2)*eye(2));

end

% Loop through layers
for il=1:length(L)

A1=eye(2)+inv(V{il+1})*V{il};
A2=eye(2)+inv(V{il+1})*V{il+2};
B1=eye(2)-inv(V{il+1})*V{il};
B2=eye(2)-inv(V{il+1})*V{il+2};


X=exp(1i*k_z(il+1)*k_0*L(il));

S{1,1}=inv(A1-X*B2*inv(A2)*X*B1)*(X*B2*inv(A2)*X*A1-B1);
S{1,2}=inv(A1-X*B2*inv(A2)*X*B1)*X*(A2-B2*inv(A2)*B2);
S{2,1}=inv(A2-X*B1*inv(A1)*X*B2)*X*(A1-B1*inv(A1)*B1);
S{2,2}=inv(A2-X*B1*inv(A1)*X*B2)*(X*B1*inv(A1)*X*A2-B2);

S_global = RH_star(S_global,S);

end


%%

% Starting parameters
k_inc = k_0*sqrt(eps_r1*mu_r1)*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]';
normal=[0;0;-1];

TE_direction=cross(k_inc,normal)/norm(cross(k_inc,normal));
nancheck=all(isnan(TE_direction));
TE_direction(isnan(TE_direction))=0;
a_TE=TE_direction+nancheck*[0,1,0]';
a_TM=cross(a_TE,k_inc)/norm(cross(a_TE,k_inc));

p=pte*a_TE+ptm*a_TM;

p=p/norm(p);

c_inc=p(1:2);
c_ref=S_global{1,1}*c_inc;
c_trn=S_global{2,1}*c_inc;

E_ref=c_ref;
E_trn=c_trn;
E_inc=p;

E_ref(3)=-(k_x*E_ref(1)+k_y*E_ref(2))/k_z(1);
E_trn(3)=-(k_x*E_trn(1)+k_y*E_trn(2))/k_z(end);
E_inc(3)=-(k_x*E_inc(1)+k_y*E_inc(2))/k_z(1);

R(iter)=dot(E_ref,E_ref)/dot(E_inc,E_inc); %maybe conj?

T(iter)=dot(E_trn,E_trn)/dot(E_inc,E_inc)*real(mu_r1/mu_r2*k_z(end)/k_z(1));

end

Abs(:,1)=1-R-T;
Abs(:,2)=R;
Abs(:,3)=T;
pll=[300:1000';300:1000';300:1000']';
area(pll,Abs)
legend("GaAs","R","T")


function eps=getdata(x,materials,name,lab)

index = find(strcmp(materials, name));
ind=find(x{1,index}(:,1)==lab);

eps=(x{1,index}(ind,2)+1i*x{1,index}(ind,3))^2;

end