function layer = apply_convolution(layer, param, iter)

for i = 1:numel(layer)
    % Start fresh at every wavelength
    layer(i).geometry.eps = layer(i).geometry.eps_struc;
    % Loop though materials in multilayer
    for j = 1:numel(layer(i).permittivities)
        % Assign permittivity at specific wavelength using interpolant
        eps = layer(i).permittivities{j}(param.wavelengthArray(iter));
        layer(i).geometry.eps(layer(i).geometry.eps_struc == j) = eps;
    end

    if iscell(layer(i).input)
        % Convolution
        layer(i).geometry.eps = conv_mat(layer(i).geometry.eps, param.P, param.Q);
        layer(i).geometry.mu_struc = eye(param.P * param.Q);

        % Generate mu and truncate mu and eps
        layer(i).geometry.mu = layer(i).geometry.mu_struc(param.tr_ind, param.tr_ind, :);
        layer(i).geometry.eps = layer(i).geometry.eps(param.tr_ind, param.tr_ind, :);
    end
end
end


function C = conv_mat(A, P, Q)
% Create the convolution matrix

p_r = -floor((P - 1) / 2) : floor((P - 1) / 2);
q_r = -floor((Q - 1) / 2) : floor((Q - 1) / 2);

[Nx, Ny, zdim] = size(A);

% Resize FFT if resolution is lower than 2x harmonics
FFT_resize_x = max(2 * P, Nx);
FFT_resize_y = max(2 * Q, Ny);

C = zeros(P * Q, P * Q, zdim);

% Middle of array
c0 = floor(P * Q / 2) + 1;
p0 = floor(FFT_resize_x / 2) + 1;
q0 = floor(FFT_resize_y / 2) + 1;

for v = 1:zdim

    B = fftshift(fft2(A(:, :, v), FFT_resize_x, FFT_resize_y)) / (Nx * Ny);

    for p_acc = p_r %p'
        for q_acc = q_r %q'
            for p_ast = p_r %p
                for q_ast = q_r %q
                    C(c0 + q_acc * P + p_acc, c0 + q_ast * P + p_ast, v) = ...
                        B(p0 + p_acc - p_ast, q0 + q_acc - q_ast);
                end
            end
        end
    end

end
end

