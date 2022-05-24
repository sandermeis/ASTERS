function [V,W,lam] = calc_VW(mu_r,eps_r,Kx,Ky)

P = calc_PQ(mu_r,eps_r,Kx,Ky);
Q = calc_PQ(eps_r,mu_r,Kx,Ky);

[W, eigval]=eig(P*Q);
lam=sqrt(eigval);
V=Q*W*inv(lam);

end

