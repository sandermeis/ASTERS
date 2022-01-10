function [eps_r,mu_r]=build_perm(n,shape,lenx,leny,P,Q,eps_lab,iter)
RoughDim=10;
Roughness=3;
r=0;

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS","GaP"];

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
    elseif shape_f(i)=="GratingX"
        eps_r{i}=conv_mat(epsmatrix(1,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)=="GratingY"
        eps_r{i}=conv_mat(epsmatrix(2,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)=="GratingXY"
        eps_r{i}=conv_mat(epsmatrix(3,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)=="Triangle"
        eps_r{i}=conv_mat(epsmatrix(4,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)=="Circle"
        eps_r{i}=conv_mat(epsmatrix(5,eps,lenx,leny),P,Q);
        mu_r{i}=mu*eye(num_H);
    elseif shape_f(i)=="Rough"
        epsprev=eps_lab{find(strcmp(materials, n_f(i-1)))}(iter);
        r=i;
        eps_r{i}=conv_mat(roughsurf(eps,epsprev,lenx,leny,RoughDim,Roughness),P,Q);
        mu_r{i}=mu*eye(num_H);
    end
    
end

if r
    eps_r=[eps_r(1:r-1),squeeze(num2cell(eps_r{r},[1 2]))',eps_r(r+1:end)];
    mu_r=[mu_r(1:r-1),repmat({mu_r{r}},1,RoughDim),mu_r(r+1:end)];
end

end