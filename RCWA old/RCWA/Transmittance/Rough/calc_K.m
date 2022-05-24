function [Kx,Ky,Kz,beta] = calc_K(k_0,theta,phi,eps_r,mu_r,P,Q,labda_x,labda_y)

%%% NEED TO USE SPARSE

M = -(P-1)/2:(P-1)/2;
N = -(Q-1)/2:(Q-1)/2;

n_inc=sqrt(eps_r{1}(1,1)*mu_r{1}(1,1));
beta(1)=k_0*n_inc*cos(phi)*sin(theta);
beta(2)=k_0*n_inc*sin(phi)*sin(theta);
beta(3)=k_0*n_inc*cos(theta);

k_x_inc=(beta(1)-2*pi*M/labda_x)/k_0;
k_y_inc=(beta(2)-2*pi*N/labda_y)/k_0;

[ky,kx] = meshgrid(k_y_inc,k_x_inc);

Kx = diag(kx(:));
Ky = diag(ky(:));
% len=length(k_x_inc);
% Kx=kron(diag(k_x_inc(:)),eye(len));
% Ky=kron(eye(len),diag(k_y_inc(:)));

for i=1:size(eps_r,2)
Kz{i}=conj(sqrt(conj(eps_r{i})*conj(mu_r{i})*eye(size(Kx))-Kx^2-Ky^2));
end

end