function [Sz, fields, Kx, Ky, Kz] = run_RCWA(layer, param, mode, progressTick)

if mode=="S"
    % Initialize results
    Sz = zeros(param.num_H, ~param.calcAllRough * numel(layer) + param.calcAllRough * numel([layer.L]) + 2, numel(param.wavelengthArray));
    Kx = zeros(param.num_H, numel(param.wavelengthArray));
    Ky = zeros(param.num_H, numel(param.wavelengthArray));
    Kz = cell(1, numel(param.wavelengthArray));
    fields = [];
    
    % Run sim
    for iWavelength = 1:numel(param.wavelengthArray)
        [Sz(:, :, iWavelength), Kx(:,iWavelength), Ky(:,iWavelength), Kz{iWavelength}] = RCWA_wl_SM(layer, param, iWavelength);
        progressTick();
    end
else
    % Initialize results
    Sz = zeros(param.num_H, ~param.calcAllRough * numel(layer) + param.calcAllRough * numel([layer.L]) + 2, numel(param.wavelengthArray));
    Kx = zeros(param.num_H, numel(param.wavelengthArray));
    Ky = zeros(param.num_H, numel(param.wavelengthArray));
    Kz = cell(1, numel(param.wavelengthArray));
    fields = cell(1, numel(param.wavelengthArray));
    field_resolution = 10;

    % Run sim
    for iWavelength = 1:numel(param.wavelengthArray)
        [Sz(:, :, iWavelength), Kx(:,iWavelength), Ky(:,iWavelength), Kz{iWavelength}, layer, param, F, X, c_im] = RCWA_wl_ETM(layer, param, iWavelength);
    
        if param.calcFields
            fields{iWavelength} = calc_fields(layer, param, iWavelength, field_resolution, F, X, c_im);
        end
        
        progressTick();
    end
end


end


function [Sz, Kx, Ky, Kz, layer, param, F, X, c_im] = RCWA_wl_ETM(layer, param, iWavelength)

%% Parameters %%

% Calculate incoming waves, convolution matrices, K matrices
param = get_sinc(param, iWavelength);
layer = apply_convolution(layer, param, iWavelength);
[layer, param] = calc_K(layer, param, iWavelength);

mu_ref = param.mu_ref(iWavelength);
mu_trn = param.mu_trn(iWavelength);

% Boundary conditions reflection region

A1 = -1i / mu_ref * (param.Kz_ref \ param.Kx) * param.Ky;
A2 = -1i / mu_ref * (param.Kz_ref \ (param.Ky * param.Ky + param.Kz_ref * param.Kz_ref));
A3 = 1i / mu_ref * (param.Kz_ref \ (param.Kx * param.Kx + param.Kz_ref * param.Kz_ref));
A4 = 1i / mu_ref * (param.Kz_ref \ param.Kx) * param.Ky;

A = [eye(2 * param.num_H);...
    A1, A2;...
    A3, A4];

% Boundary conditions transmission region

B1 = 1i / mu_trn * (param.Kz_trn \ param.Kx) * param.Ky;
B2 = 1i / mu_trn * (param.Kz_trn \ (param.Ky * param.Ky + param.Kz_trn * param.Kz_trn));
B3 = -1i / mu_trn * (param.Kz_trn \ (param.Kx * param.Kx + param.Kz_trn * param.Kz_trn));
B4 = -1i / mu_trn * (param.Kz_trn \ param.Kx) * param.Ky;

B=[eye(2 * param.num_H);...
    B1, B2;...
    B3, B4];

% Set equal athe the boundaries
% s = incoming
% r = reflection
% t = transmission
%
% matrix composition:
% 
% 1st quarter electric field in x
% 2nd quarter electric field in y
% 3rd quarter magnetic field in x
% 4th quarter magnetic field in y
%
%  ref              int               trn
%      |                           |      
%      |                           |      
% s+Ar | F1[I, 0]c1     F1[X1,0]c1 |  Bt
%      |   [0,X1]         [0, I]   |       
%      |                           |       


%% Loop backwards through layers %%

for i = numel(layer):-1:1
    % Loop backwards through sublayers
    for j = size(layer(i).geometry.eps, 3):-1:1
        
        % Calculate V, W matrices and eigenvalue
        [V, W, lam] = calc_VW(layer(i).geometry.mu, layer(i).geometry.eps(:, :, j), param.Kx, param.Ky);
        
        % Assemble field matrix
        F{i,j} = [W, W; -V, V];

        % Assemble propagation matrix
        X{i,j} = expm(-lam * param.k_0(iWavelength) * layer(i).L(j));

        % Get a and b matrices, [a; b] = F^(-1) B
        temp1 = F{i,j} \ B;

        a{i,j} = temp1(1:2 * param.num_H, :);
        b{i,j} = temp1(2 * param.num_H + 1:end, :);

        % Effective 'B' of lower layer
        B = F{i,j} * [eye(2*param.num_H), zeros(2*param.num_H); zeros(2*param.num_H), X{i, j}] * [eye(2*param.num_H); b{i, j} * (a{i, j} \ X{i, j})];
    end
end

% Now B is 'complete', solve for reflection region boundary conditions

temp2 = [-A, B] \ param.s_inc;
r = temp2(1:2 * param.num_H);
tN{1} = temp2(2 * param.num_H + 1:end);

%% Looping forward: with t1, a and X, now all t are solved recursively %%

k = 1;
for i = 1:numel(layer)
    for j = 1:size(layer(i).geometry.eps, 3)
        temp3 = a{i,j} \ X{i,j};
        tN{k + 1} = temp3 * tN{k};
        c_im{i, j} = [eye(2 * param.num_H); b{i, j} * temp3] * tN{k};

        % calcAllRough 1: All sublayer fields are calculated
        if param.calcAllRough
            % Calculate fields
            field_begin = F{i, j} * [eye(2 * param.num_H), zeros(2 * param.num_H); zeros(2 * param.num_H), X{i, j}] * c_im{i, j};
            field_end = F{i, j} * [X{i, j}, zeros(2 * param.num_H); zeros(2 * param.num_H), eye(2 * param.num_H)] * c_im{i, j};

            % Calculate poynting vector
            Sz(:, k) = (Sz_field(field_begin) - Sz_field(field_end)) / sum(Sz_field(param.s_inc));
        end
        k = k + 1;

    end

    % calcAllRough 0: Only first and last sublayer fields are calculated
    if ~param.calcAllRough
        % Calculate field at start of first layer
        field_begin = F{i, 1} * [eye(2 * param.num_H), zeros(2 * param.num_H); zeros(2 * param.num_H), X{i, 1}] * c_im{i, 1};
        % Calculate field at end of last layer
        field_end = F{i, j} * [X{i,j}, zeros(2 * param.num_H); zeros(2 * param.num_H), eye(2 * param.num_H)] * c_im{i, j};

        % Calculate poynting vector
        Sz(:, i) = (Sz_field(field_begin) - Sz_field(field_end)) / sum(Sz_field(param.s_inc));
    end

end


%% Diffraction efficiencies and return parameters %%

t = temp3 * tN{end - 1};

r(2 * param.num_H + 1:3 * param.num_H) = - param.Kz_ref \ (param.Kx * r(1:param.num_H) + param.Ky * r(param.num_H + 1:2 * param.num_H));
t(2 * param.num_H + 1:3 * param.num_H) = - param.Kz_trn \ (param.Kx * t(1:param.num_H) + param.Ky * t(param.num_H + 1:2 * param.num_H));

% Reflection region
rmag = abs(r(1:param.num_H)).^2 + abs(r(param.num_H + 1:2 * param.num_H)).^2 + abs(r(2 * param.num_H + 1:end)).^2;
Sz(:, end + 1) = real(param.Kz_ref * rmag / param.beta(3));

% Transmission region
tmag = abs(t(1:param.num_H)).^2 + abs(t(param.num_H + 1:2 * param.num_H)).^2 + abs(t(2 * param.num_H + 1:end)).^2;
Sz(:, end + 1) = real(param.Kz_trn * tmag / param.beta(3));

% Return wave vectors
Kx = diag(param.Kx);
Ky = diag(param.Ky);

for index = 1:numel(layer)
    Kz{index} = 0;%layer(index).geometry.Kz;
end

end

% Scattering matrix formalism %% UNSTABLE FOR TEXTURED LAYERS
function [Sz, Kx, Ky, Kz] = RCWA_wl_SM(layer, param, iWavelength)

%% Parameters %%

% Calculate incoming waves, convolution matrices, K matrices
param = get_sinc(param, iWavelength);
layer = apply_convolution(layer, param, iWavelength);
[layer, param] = calc_K(layer, param, iWavelength);

mu_ref = param.mu_ref(iWavelength);
mu_trn = param.mu_trn(iWavelength);
eps_ref = param.mu_ref(iWavelength);
eps_trn = param.mu_trn(iWavelength);

% Free space
SI = {zeros(2 * param.num_H), eye(2 * param.num_H); eye(2 * param.num_H), zeros(2 * param.num_H)};
[V0, W0, ~] = calc_VW(eye(size(param.Kx)), eye(size(param.Kx)), param.Kx, param.Ky);

% Calculate V, W matrices and eigenvalue
k = 1;
[V_ref, W_ref, ~] = calc_VW(mu_ref, eps_ref, param.Kx, param.Ky);

A = inv(W_ref) * W0 + inv(V_ref) * V0;
B = inv(W_ref) * W0 - inv(V_ref) * V0;

S{1, 1} = inv(A - B * inv(A) * B) * (B * inv(A) * A - B);
S{2, 2} = S{1, 1};
S{1, 2} = inv(A - B * inv(A) * B) * (A - B * inv(A) * B);
S{2, 1} = S{1, 2};

Sg{k} = RH_star(SI, S);

for i = 1:numel(layer)
    % Loop through sublayers
    for j = 1:size(layer(i).geometry.eps, 3)
        k = k + 1;
        [V{i, j}, W{i, j}, lam] = calc_VW(layer(i).geometry.mu, layer(i).geometry.eps(:, :, j), param.Kx, param.Ky);


        A = inv(W{i, j}) * W0 + inv(V{i, j}) * V0;
        B = inv(W{i, j}) * W0 - inv(V{i, j}) * V0;
        
        % Assemble propagation matrix
        X{i, j} = expm(-lam * param.k_0(iWavelength) * layer(i).L(j));
        
        S{1, 1} = inv(A - X{i, j} * B * inv(A) * X{i, j} * B) * (X{i, j} * B * inv(A) * X{i, j} * A - B);
        S{2, 2} = S{1, 1};
        S{1, 2} = inv(A - X{i, j} * B * inv(A) * X{i, j} * B) * X{i, j} * (A - B * inv(A) * B);
        S{2, 1} = S{1, 2};
        
        Sg{k} = RH_star(Sg{k - 1}, S);
    end
end

k = k + 1;
[V_trn, W_trn, ~] = calc_VW(mu_trn, eps_trn, param.Kx, param.Ky);

A = inv(W_trn) * W0 + inv(V_trn) * V0;
B = inv(W_trn) * W0 - inv(V_trn) * V0;

S{1, 1} = inv(A - B * inv(A) * B) * (B * inv(A) * A - B);
S{2, 2} = S{1, 1};
S{1, 2} = inv(A - B * inv(A) * B) * (A - B * inv(A) * B);
S{2, 1} = S{1, 2};

Sg{k} = RH_star(Sg{k-1}, S);

c_inc = inv(W_ref) * param.s_inc(1:2 * param.num_H);
c_ref = Sg{end}{1, 1} * c_inc;
c_trn = Sg{end}{2, 1} * c_inc;

% Fields
for i = numel(layer):-1:1
    for j = size(layer(i).geometry.eps, 3):-1:1

        k = k - 1;
        c_mn = inv(Sg{k - 1}{1, 2}) * (c_ref - Sg{k - 1}{1, 1} * c_inc);
        c_pl = Sg{k - 1}{2, 1} * c_inc + Sg{k - 1}{2, 2} * c_mn;
        A = inv(W{i, j}) * W0 + inv(V{i, j}) * V0;
        B = inv(W{i, j}) * W0 - inv(V{i, j}) * V0;
        c{i, j} = [0.5*A, 0.5*B; 0.5*B, 0.5*A] * [c_pl; c_mn];

        % calcAllRough 1: All sublayer fields are calculated
        if param.calcAllRough
            % Calculate fields
            field_begin = [W{i, j}, W{i, j}; -V{i, j}, V{i, j}] * c{i, j};
            field_end = [W{i, j}, W{i, j}; -V{i, j}, V{i, j}] * [X{i, j}, zeros(2 * param.num_H); zeros(2 * param.num_H), inv(X{i, j})] * c{i,j};
            % Calculate poynting vector
            Sz(:, k-1) = (Sz_field(field_begin) - Sz_field(field_end)) / sum(Sz_field(param.s_inc));
        end
        

    end

    % calcAllRough 0: Only first and last sublayer fields are calculated
    if ~param.calcAllRough
        % Calculate field at start of first layer
        field_begin = [W{i, 1}, W{i, 1}; -V{i, 1}, V{i, 1}] * c{i, 1};
        % Calculate field at end of last layer
        field_end = [W{i, j}, W{i, j}; -V{i, j}, V{i, j}] * [X{i, j}, zeros(2 * param.num_H); zeros(2 * param.num_H), inv(X{i, j})] * c{i,j};

        % Calculate poynting vector
        Sz(:, i) = (Sz_field(field_begin) - Sz_field(field_end)) / sum(Sz_field(param.s_inc));
    end

end


%% Diffraction efficiencies and return parameters %%

r = W_ref * c_ref;
t = W_trn * c_trn;

r(2 * param.num_H + 1:3 * param.num_H) = - param.Kz_ref \ (param.Kx * r(1:param.num_H) + param.Ky * r(param.num_H + 1:2 * param.num_H));
t(2 * param.num_H + 1:3 * param.num_H) = - param.Kz_trn \ (param.Kx * t(1:param.num_H) + param.Ky * t(param.num_H + 1:2 * param.num_H));

% Reflection region
rmag = abs(r(1:param.num_H)).^2 + abs(r(param.num_H + 1:2 * param.num_H)).^2 + abs(r(2 * param.num_H + 1:end)).^2;
Sz(:, end + 1) = real(param.Kz_ref * rmag / param.beta(3));

% Transmission region
tmag = abs(t(1:param.num_H)).^2 + abs(t(param.num_H + 1:2 * param.num_H)).^2 + abs(t(2 * param.num_H + 1:end)).^2;
Sz(:, end + 1) = real(param.Kz_trn * tmag / param.beta(3));

% Return wave vectors
Kx = diag(param.Kx);
Ky = diag(param.Ky);

for index = 1:numel(layer)
    Kz{index} = 0;%layer(index).geometry.Kz;
end
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
    PQ = 1 / mu * [Kx * Ky, eps * mu * eye(size(Kx)) - Kx * Kx;
        Ky * Ky - eps * mu * eye(size(Kx)), -Ky * Kx];
else
    PQ = [Kx * inv(mu) * Ky, eps - Kx * inv(mu) * Kx;
        Ky * inv(mu) * Ky - eps, - Ky * inv(mu) * Kx];
end

end


function sz = Sz_field(f)
F = reshape(f, [], 4);
sz = 0.5 * real(F(:, 1) .* conj(1i * F(:, 4)) - F(:, 2) .* conj(1i * F(:, 3)));
end


function fields = calc_fields(layer, param, iWavelength, field_resolution, F, X, c_im)

fields = cell(1, numel(layer));

% Add phase of incoming field
[x_grid, y_grid] = ndgrid(linspace(0, param.size, param.res));
phase = exp(1i * param.k_0(iWavelength) * (param.beta(1) * x_grid + param.beta(2) * y_grid));

% Initialize field arrays
%Ex = zeros(param.res, param.res, sum([layer.L]));
%Ey = Ex; Ez = Ex; Hx = Ex; Hy = Ex; Hz = Ex;

% Loop through layers
total_index = 1;

for i = 1:numel(layer)
    num_zdim = size(layer(i).geometry.eps, 3);
    % Loop through sublayers
    for j = 1:num_zdim

        % Local grids
        d = layer(i).L(j);
        field_resolution = d + 1;
        z_grid = linspace(0, d, field_resolution);

        for z = 1:numel(z_grid) - 1
            
            f = F{i, j} * [diag(diag(X{i, j}).^(z_grid(z) / d)), zeros(2 * param.num_H);...
                        zeros(2 * param.num_H), diag(diag(X{i, j}).^(-(z_grid(z) - d) / d))] * c_im{i, j};

            f2 = reshape(f,[],4);
            ezf = -1i * inv(layer(i).geometry.eps(:, :, j)) * (param.Kx * f2(:,4) - param.Ky * f2(:, 3));
            hzf = -1i * inv(layer(i).geometry.mu) * (param.Kx * f2(:, 2) - param.Ky * f2(:, 1));

            Ex(:, :, total_index) = phase .* ifftn(reshape(f2(:, 1), param.Hmax, param.Hmax), [param.res, param.res]);
            Ey(:, :, total_index) = phase .* ifftn(reshape(f2(:, 2), param.Hmax, param.Hmax), [param.res, param.res]);
            Ez(:, :, total_index) = phase .* ifftn(reshape(ezf, param.Hmax, param.Hmax), [param.res, param.res]);

            Hx(:, :, total_index) = phase .* ifftn(reshape(f2(:, 3), param.Hmax, param.Hmax), [param.res, param.res]);
            Hy(:, :, total_index) = phase .* ifftn(reshape(f2(:, 4), param.Hmax, param.Hmax), [param.res, param.res]);


            Hz(:, :, total_index) = phase .* ifftn(reshape(hzf, param.Hmax, param.Hmax), [param.res, param.res]);
            total_index = total_index + 1;
        end

    end
    
end

fields = {Ex, Ey, Ez};%,Hx,Hy,Hz};

end


function S_AB = RH_star(S_A, S_B)
% Redheffer Star product
% Input is cell array of matrices.
sz = size(S_A{1, 1});
D = S_A{1, 2} / (eye(sz) - S_B{1, 1} * S_A{2, 2});
F = S_B{2, 1} / (eye(sz) - S_A{2, 2} * S_B{1, 1});

S_AB{1, 1} = S_A{1, 1} + D * S_B{1, 1} * S_A{2, 1};
S_AB{1, 2} = D * S_B{1, 2};
S_AB{2, 1} = F * S_A{2, 1};
S_AB{2, 2} = S_B{2, 2} + F * S_A{2, 2} * S_B{1, 2};
end