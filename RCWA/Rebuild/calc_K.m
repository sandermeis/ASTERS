function [layer, param] = calc_K(layer, param, iter)

M = -(param.P-1)/2:(param.P-1)/2;
N = -(param.Q-1)/2:(param.Q-1)/2;

k_x_inc = param.beta(1)-2*pi*M/(param.size_x*param.k_0(iter));
k_y_inc = param.beta(2)-2*pi*N/(param.size_y*param.k_0(iter));

[ky, kx] = meshgrid(k_y_inc, k_x_inc);

param.Kx = diag(kx(:));
param.Ky = diag(ky(:));

param.Kx = param.Kx(param.tr_ind,param.tr_ind);
param.Ky = param.Ky(param.tr_ind,param.tr_ind);

% len=length(k_x_inc);
% Kx=kron(diag(k_x_inc(:)),eye(len));
% Ky=kron(eye(len),diag(k_y_inc(:)));

eps_ref = param.eps_ref(iter);
eps_trn = param.eps_trn(iter);
mu_ref = param.mu_ref(iter);
mu_trn = param.mu_trn(iter);

param.Kz_ref = conj(sqrt(conj(eps_ref)*conj(mu_ref)*eye(size(param.Kx))-param.Kx^2-param.Ky^2));
param.Kz_trn = conj(sqrt(conj(eps_trn)*conj(mu_trn)*eye(size(param.Kx))-param.Kx^2-param.Ky^2));

for i = 1:numel(layer)
    layer(i).geometry.Kz = zeros(size(param.Kx));
    for j = 1:size(layer(i).geometry.eps,3)
        layer(i).geometry.Kz(:,:,j) = conj(sqrt(conj(layer(i).geometry.eps(:,:,j))*conj(layer(i).geometry.mu(:,:))*eye(size(param.Kx))-param.Kx^2-param.Ky^2));
    end
end

end