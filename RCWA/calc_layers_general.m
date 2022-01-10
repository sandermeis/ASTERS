function [Sg,W,V,lam] = calc_layers_general(num_H,mu_r,eps_r,Kx,Ky,Lk)

Sg{1} = {zeros(2*num_H),eye(2*num_H);eye(2*num_H),zeros(2*num_H)};

for i=1:length(eps_r)
[V{i},W{i},lam{i}] = calc_VW(mu_r{i},eps_r{i},Kx,Ky);
end

for i=2:(length(eps_r)-1)

wi=inv(W{i});
vi=inv(V{i});
A1=wi*W{i-1}+vi*V{i-1};
A2=wi*W{i+1}+vi*V{i+1};
B1=wi*W{i-1}-vi*V{i-1};
B2=wi*W{i+1}-vi*V{i+1};

X=expm(-lam{i}*Lk(i));

S{1,1}=inv(A1-X*B2*inv(A2)*X*B1)*(X*B2*inv(A2)*X*A1-B1);
S{1,2}=inv(A1-X*B2*inv(A2)*X*B1)*X*(A2-B2*inv(A2)*B2);
S{2,1}=inv(A2-X*B1*inv(A1)*X*B2)*X*(A1-B1*inv(A1)*B1);
S{2,2}=inv(A2-X*B1*inv(A1)*X*B2)*(X*B1*inv(A1)*X*A2-B2);

Sg{i} = RH_star(Sg{i-1},S);
end

end

