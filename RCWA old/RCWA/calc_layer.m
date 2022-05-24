function [S,W,V,lam] = calc_layer(V_0,W_0,mu_r,eps_r,Kx,Ky,Lk)

P = calc_PQ(mu_r,eps_r,Kx,Ky);
Q = calc_PQ(eps_r,mu_r,Kx,Ky);

[W, eigval]=eig(P*Q);
lam=sqrt(eigval);
V=Q*W*inv(lam);

A=inv(W)*W_0+inv(V)*V_0;
B=inv(W)*W_0-inv(V)*V_0;

X = expm(-lam*Lk);

S{1,1}=inv(A-X*B*inv(A)*X*B)*(X*B*inv(A)*X*A-B);
S{2,2}=S{1,1};

S{1,2}=inv(A-X*B*inv(A)*X*B)*X*(A-B*inv(A)*B);
S{2,1}=S{1,2};

end

