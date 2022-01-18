function [r,t,J] = calc_layers_transmittance(num_H,mu_r0,eps_r0,Kx,Ky,Kz,k_0,s_inc,L,shape,roughdim)

indicator=ones(size(L));
for i=find(shape==6|shape==7)
    L(1:roughdim(i),i)=L(1,i)/roughdim(i)*ones(1,roughdim(i));
    indicator(1:roughdim(i),i)=2*ones(1,roughdim(i));
end
indicator=indicator(indicator~=0).';
indicator=indicator-1;
L=L(L~=0).';

mu_r_ref=mu_r0{1};
mu_r_trn=mu_r0{end};

eps_r=eps_r0(2:end-1);
mu_r=mu_r0(2:end-1);

N=length(eps_r);

A1=-1i/mu_r_ref*inv(Kz{1})*Kx*Ky;
A2=-1i/mu_r_ref*inv(Kz{1})*(Ky*Ky+Kz{1}*Kz{1});
A3=1i/mu_r_ref*inv(Kz{1})*(Kx*Kx+Kz{1}*Kz{1});
A4=1i/mu_r_ref*inv(Kz{1})*Kx*Ky;

B1=1i/mu_r_trn*inv(Kz{end})*Kx*Ky;
B2=1i/mu_r_trn*inv(Kz{end})*(Ky*Ky+Kz{end}*Kz{end});
B3=-1i/mu_r_trn*inv(Kz{end})*(Kx*Kx+Kz{end}*Kz{end});
B4=-1i/mu_r_trn*inv(Kz{end})*Kx*Ky;

A=[eye(2*num_H);...
    A1,A2;...
    A3,A4];
B=[eye(2*num_H);...
    B1,B2;...
    B3,B4];

for i=N:-1:1
[V,W,lam] = calc_VW(mu_r{i},eps_r{i},Kx,Ky);

F{i} = [W,W;-V,V];
X{i}=expm(-lam*k_0*L(i));

zz = mat2cell(inv(F{i})*B,2*num_H*ones(1,2));

a{i}=zz{1};
b{i}=zz{2};


B = F{i}*[eye(2*num_H),zeros(2*num_H);zeros(2*num_H),X{i}]*[eye(2*num_H);b{i}*inv(a{i})*X{i}];
end

zz=inv([-A, B]);
res=mat2cell(zz*s_inc,[2*num_H,2*num_H]);
r=mat2cell(res{1},[num_H,num_H]);

tN{1}=res{2};

j=1;

% for i=1:N
%     tN{i+1}=inv(a{i})*X{i}*tN{i};
%     c_im=[eye(2*num_H);b{i}*inv(a{i})*X{i}]*tN{i};
%     
%     field_begin=F{i}*[eye(2*num_H),zeros(2*num_H);zeros(2*num_H),X{i}]*c_im;
%     field_end=F{i}*[X{i},zeros(2*num_H);zeros(2*num_H),eye(2*num_H)]*c_im;
%     J{i}=(sum(Sz_field(field_begin))-sum(Sz_field(field_end)))/sum(Sz_field(s_inc));
%     
% end

for i=1:N
    
    tN{i+1}=inv(a{i})*X{i}*tN{i};
    c_im=[eye(2*num_H);b{i}*inv(a{i})*X{i}]*tN{i};
    
    if indicator(i)==1
        try
            if indicator(i-1)==0
                field_begin=F{i}*[eye(2*num_H),zeros(2*num_H);zeros(2*num_H),X{i}]*c_im;
            elseif indicator(i+1)==0
                field_end=F{i}*[X{i},zeros(2*num_H);zeros(2*num_H),eye(2*num_H)]*c_im;
                J{j}=(sum(Sz_field(field_begin))-sum(Sz_field(field_end)))/sum(Sz_field(s_inc));
                j=j+1;
            end
        catch ME
            if strcmp(ME.message,'Array indices must be positive integers or logical values.')
                field_begin=F{i}*[eye(2*num_H),zeros(2*num_H);zeros(2*num_H),X{i}]*c_im;
            elseif strcmp(ME.identifier,'MATLAB:badsubscript')
                field_end=F{i}*[X{i},zeros(2*num_H);zeros(2*num_H),eye(2*num_H)]*c_im;
                J{j}=(sum(Sz_field(field_begin))-sum(Sz_field(field_end)))/sum(Sz_field(s_inc));
                j=j+1;
            else
                rethrow(ME)
            end
        end
    else
        field_begin=F{i}*[eye(2*num_H),zeros(2*num_H);zeros(2*num_H),X{i}]*c_im;
        field_end=F{i}*[X{i},zeros(2*num_H);zeros(2*num_H),eye(2*num_H)]*c_im;
        J{j}=(sum(Sz_field(field_begin))-sum(Sz_field(field_end)))/sum(Sz_field(s_inc));
        j=j+1;
    end
    
end

t=mat2cell(inv(a{N})*X{N}*tN{N},[num_H,num_H]);

end

