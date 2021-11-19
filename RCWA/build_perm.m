function c=build_perm(ref_medium,trn_medium,n,shape,eps_lab,lab)

n_f=[ref_medium,n,trn_medium];
shape_f=["Uniform",shape,"Uniform"];
len_n=length(n_f);

for i=1:len_n
index = find(strcmp(materials, n_f(i))); %find tab index
n_eps(i)=eps_lab{index}(lab);
end

% eps_input=epsmatrix(getdata(x,materials,n(1),lam0),lenx,leny);
% eps_material=conv_mat(eps_input,P,Q);
% mu_material=1;
% 
% eps_trn=1;
% mu_trn=1;
% eps_r   = {eps_ref*eye(num_H), eps_material, eps_trn*eye(num_H)};
% mu_r    = {mu_ref*eye(num_H), mu_material*eye(num_H), mu_trn*eye(num_H)};
% L       = [0, 300, 0];

end
