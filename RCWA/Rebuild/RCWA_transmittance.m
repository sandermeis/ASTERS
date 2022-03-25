function Sz = RCWA_transmittance(layer, param)

Sz = zeros(param.num_H,~param.calcAllRough*numel(layer)+param.calcAllRough*numel([layer.L])+2,numel(param.wavelengthArray));

for iWavelength = 1:numel(param.wavelengthArray)
    Sz(:, :, iWavelength) = RCWA_wl(layer, param, iWavelength);
end

end


function Sz = RCWA_wl(layer, param, iWavelength)

param = get_sinc(param, iWavelength);
layer = apply_convolution(layer, param, iWavelength);
[layer, param] = calc_K(layer, param, iWavelength);

mu_ref = param.mu_ref(iWavelength);
mu_trn = param.mu_trn(iWavelength);

A1 = -1i/mu_ref*(param.Kz_ref\param.Kx)*param.Ky;
A2 = -1i/mu_ref*(param.Kz_ref\(param.Ky*param.Ky+param.Kz_ref*param.Kz_ref));
A3 = 1i/mu_ref*(param.Kz_ref\(param.Kx*param.Kx+param.Kz_ref*param.Kz_ref));
A4 = 1i/mu_ref*(param.Kz_ref\param.Kx)*param.Ky;

B1 = 1i/mu_trn*(param.Kz_trn\param.Kx)*param.Ky;
B2 = 1i/mu_trn*(param.Kz_trn\(param.Ky*param.Ky+param.Kz_trn*param.Kz_trn));
B3 = -1i/mu_trn*(param.Kz_trn\(param.Kx*param.Kx+param.Kz_trn*param.Kz_trn));
B4 = -1i/mu_trn*(param.Kz_trn\param.Kx)*param.Ky;

A=[eye(2*param.num_H);...
    A1,A2;...
    A3,A4];
B=[eye(2*param.num_H);...
    B1,B2;...
    B3,B4];

for i = numel(layer):-1:1
    for j = size(layer(i).geometry.eps,3):-1:1

        [V, W, lam] = calc_VW(layer(i).geometry.mu(:,:),layer(i).geometry.eps(:,:,j),param.Kx,param.Ky);
        F{i,j} = [W, W; -V, V];
        X{i,j} = expm(-lam * param.k_0(iWavelength) * layer(i).L(j));

        temp = F{i,j}\B;

        a{i,j} = temp(1:2*param.num_H,:);
        b{i,j} = temp(2*param.num_H+1:end,:);

        B = F{i,j} * [eye(2*param.num_H),zeros(2*param.num_H);zeros(2*param.num_H),X{i,j}] * [eye(2*param.num_H);b{i,j}*(a{i,j}\X{i,j})];
    end
end

temp = [-A, B]\param.s_inc;
r = temp(1:2*param.num_H);
tN{1} = temp(2*param.num_H+1:end);

k = 1;
if param.calcAllRough
    Sz = zeros(param.num_H,numel([layer.L]));
else
    Sz = zeros(param.num_H,numel(layer));
end

for i = 1:numel(layer)
    for j = 1:size(layer(i).geometry.eps,3)
        temp = a{i,j}\X{i,j};
        tN{k+1} = temp*tN{k};
        c_im = [eye(2*param.num_H);b{i,j}*temp]*tN{k};

        if param.calcAllRough
            field_begin = F{i,j}*[eye(2*param.num_H),zeros(2*param.num_H);zeros(2*param.num_H),X{i,j}]*c_im;
            field_end = F{i,j}*[X{i,j},zeros(2*param.num_H);zeros(2*param.num_H),eye(2*param.num_H)]*c_im;
            Sz(:,k) = (Sz_field(field_begin)-Sz_field(field_end))/sum(Sz_field(param.s_inc));
        elseif j==1
            field_begin = F{i,j}*[eye(2*param.num_H),zeros(2*param.num_H);zeros(2*param.num_H),X{i,j}]*c_im;
        end
        k=k+1;
    end
    if ~param.calcAllRough
        field_end = F{i,j}*[X{i,j},zeros(2*param.num_H);zeros(2*param.num_H),eye(2*param.num_H)]*c_im;
        Sz(:,i) = (Sz_field(field_begin)-Sz_field(field_end))/sum(Sz_field(param.s_inc));
    end
end

t = temp*tN{end-1};

r(2*param.num_H+1:3*param.num_H) = -param.Kz_ref\(param.Kx*r(1:param.num_H)+param.Ky*r(param.num_H+1:2*param.num_H));
t(2*param.num_H+1:3*param.num_H) = -param.Kz_trn\(param.Kx*t(1:param.num_H)+param.Ky*t(param.num_H+1:2*param.num_H));

rmag = abs(r(1:param.num_H)).^2+abs(r(param.num_H+1:2*param.num_H)).^2+abs(r(2*param.num_H+1:end)).^2;
Sz(:,end+1) = real(param.Kz_ref*rmag/param.beta(3));

tmag = abs(t(1:param.num_H)).^2+abs(t(param.num_H+1:2*param.num_H)).^2+abs(t(2*param.num_H+1:end)).^2;
Sz(:,end+1) = real(param.Kz_trn*tmag/param.beta(3));
end


function [V, W, lam] = calc_VW(mu_r, eps_r, Kx, Ky)

P = calc_PQ(mu_r, eps_r, Kx, Ky);
Q = calc_PQ(eps_r, mu_r, Kx, Ky);

[W, eigval] = eig(P * Q);
lam = sqrt(eigval);
V = Q * W * inv(lam);

end


function PQ = calc_PQ(eps, mu, Kx, Ky)

%%%                   P
%%%
%%%| Kx eps^-1 Ky     | mu - Kx eps^-1 Kx |
%%%| Ky eps^-1 Ky - mu|    - Ky eps^-1 Kx |

%%%                   Q
%%%
%%%| Kx mu^-1 Ky      | eps - Kx mu^-1 Kx |
%%%| Ky mu^-1 Ky - eps|     - Ky mu^-1 Kx |

if and(isscalar(eps),isscalar(mu))
    PQ = 1/mu*[Kx*Ky, eps*mu*eye(size(Kx))-Kx*Kx;
        Ky*Ky-eps*mu*eye(size(Kx)), -Ky*Kx];
else
    PQ = [Kx*inv(mu)*Ky, eps-Kx*inv(mu)*Kx;
        Ky*inv(mu)*Ky-eps, -Ky*inv(mu)*Kx];
end

end


function sz = Sz_field(f)
F = reshape(f,[],4);
sz = 0.5*real(F(:,1).*conj(1i*F(:,4))-F(:,2).*conj(1i*F(:,3)));
end