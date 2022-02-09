function [eps_r,mu_r]=build_perm(n,shape,lenx,leny,P,Q,eps_lab,iter,roughdim,roughness)
r=0;

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS","GaP","AlInP"];

ref_medium="Air";
trn_medium="Air";

num_H=P*Q;

n_f=[ref_medium,n,trn_medium];
shape_f=[0,shape,0];
roughdim_f=[0,roughdim,0];
len_n=length(n_f);

for i=1:len_n
    index = find(strcmp(materials, n_f(i))); %find tab index
    eps=eps_lab{index}(iter);
    mu=1;
    
    if shape_f(i)==0
        eps_r{i}=eps;
        mu_r{i}=mu;
    elseif shape_f(i)==1
        eps_r{i}=conv_mat(epsmatrix(1,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)==2
        eps_r{i}=conv_mat(epsmatrix(2,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)==3
        eps_r{i}=conv_mat(epsmatrix(3,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)==4
        eps_r{i}=conv_mat(epsmatrix(4,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)==5
        eps_r{i}=conv_mat(epsmatrix(5,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)==6
        epsprev=eps_lab{find(strcmp(materials, n_f(i-1)))}(iter);
        eps_r{i}=conv_mat(roughsurf(eps,epsprev,lenx,leny,roughdim_f(i),roughness),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)==7
        epsprev=eps_lab{find(strcmp(materials, n_f(i-1)))}(iter);
        eps_r{i}=conv_mat(realsurf(eps,epsprev,roughdim_f(i)),P,Q);
        mu_r{i}=mu*eye(num_H);
    end
    
end

if any(find(shape_f==6|shape_f==7))
for i=find(shape_f==6|shape_f==7)
    eps_r(1:roughdim_f(i),i)=squeeze(num2cell(eps_r{i},[1 2]))';
    mu_r(1:roughdim_f(i),i)=repmat({mu_r{i}},1,roughdim_f(i));
end
eps_r=eps_r(~cellfun(@isempty,eps_r)).';
mu_r=mu_r(~cellfun(@isempty,mu_r)).';
end

end