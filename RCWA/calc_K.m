function [Kx,Ky,Kz] = calc_K(lam0,theta,phi,eps_r,mu_r,M,N,labda_x,labda_y)

k_0=2*pi/lam0;
n_inc=sqrt(eps_r{1}(1,1)*mu_r{1}(1,1));
k_x_inc=n_inc*sin(theta)*cos(phi)-2*pi*M/(k_0*labda_x);
k_y_inc=n_inc*sin(theta)*sin(phi)-2*pi*N/(k_0*labda_y);

len=length(k_x_inc);
Kx=sparse(kron(diag(k_x_inc(:)),eye(len)));
Ky=sparse(kron(eye(len),diag(k_y_inc(:))));

for i=1:size(eps_r,2)
Kz{i}=conj(sqrt(eps_r{i}*mu_r{i}-Kx^2-Ky^2));
end

end