function [eps_r,mu_r]=build_perm(n,shape,lenx,leny,P,Q,eps_lab,iter)

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS"];

ref_medium="Air";
trn_medium="Air";

num_H=P*Q;

n_f=[ref_medium,n,trn_medium];
shape_f=["Uniform",shape,"Uniform"];
len_n=length(n_f);

for i=1:len_n
    index = find(strcmp(materials, n_f(i))); %find tab index
    eps=eps_lab{index}(iter);
    mu=1;
    
    if shape_f(i)=="Uniform"
        eps_r{i}=eps;
        mu_r{i}=mu;
    elseif shape_f(i)=="Grating"
        eps_r{i}=conv_mat(epsmatrix(1,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)=="Triangle"
        eps_r{i}=conv_mat(epsmatrix(2,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)=="Circle"
        eps_r{i}=conv_mat(epsmatrix(3,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H); 
    end
    
end
