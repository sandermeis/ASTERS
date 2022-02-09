function [layer,device] = calc_K(layer,device,input_wave,iter)

device.k_0i = device.k_0(iter);

input_wave.s_inc(2*device.num_H+1:end) = device.k_0(iter)*input_wave.s_inc(2*device.num_H+1:end);

M = -(device.P-1)/2:(device.P-1)/2;
N = -(device.Q-1)/2:(device.Q-1)/2;

k_x_inc = input_wave.beta(1)-2*pi*M/(device.size_x*device.k_0(iter));
k_y_inc = input_wave.beta(2)-2*pi*N/(device.size_y*device.k_0(iter));

[ky,kx] = meshgrid(k_y_inc,k_x_inc);

device.Kx = diag(kx(:));
device.Ky = diag(ky(:));

device.Kx = device.Kx(device.tr_ind,device.tr_ind);
device.Ky = device.Ky(device.tr_ind,device.tr_ind);

% len=length(k_x_inc);
% Kx=kron(diag(k_x_inc(:)),eye(len));
% Ky=kron(eye(len),diag(k_y_inc(:)));

device.Kz_ref = conj(sqrt(conj(device.eps_ref)*conj(device.mu_ref)*eye(size(device.Kx))-device.Kx^2-device.Ky^2));
device.Kz_trn = conj(sqrt(conj(device.eps_trn)*conj(device.mu_trn)*eye(size(device.Kx))-device.Kx^2-device.Ky^2));

for i = 1:numel(layer)
    layer(i).geometry.Kz = zeros(size(device.Kx));
    for j = 1:size(layer(i).geometry.eps,3)
        layer(i).geometry.Kz(:,:,j) = conj(sqrt(conj(layer(i).geometry.eps(:,:,j))*conj(layer(i).geometry.mu(:,:))*eye(size(device.Kx))-device.Kx^2-device.Ky^2));
    end
end

end