function [Kx,Ky,Kz] = calc_K(k_0,theta,phi,eps_r,mu_r,M,N,labda_x,labda_y)

n_inc=sqrt(eps_r{1}(1,1)*mu_r{1}(1,1));
k_x_inc=n_inc*sin(theta)*cos(phi)-2*pi*M/(k_0*labda_x);
k_y_inc=n_inc*sin(theta)*sin(phi)-2*pi*N/(k_0*labda_y);

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