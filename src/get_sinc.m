function param = get_sinc(param, iter)

% Jones vector
J = [param.pTM; param.pTE; 0];

% Force unity when pTE and pTM do not sum to 1
J = J / norm(J);

% Phi is angle from x vector
% Theta is angle from surface normal

% Incoming wavevector
n_inc = sqrt(param.eps_ref(iter) * param.mu_ref(iter));
param.beta(1) = n_inc * sin(param.theta) * cos(param.phi);
param.beta(2) = n_inc * sin(param.theta) * sin(param.phi);
param.beta(3) = n_inc * cos(param.theta);

%% Vector way
% % Normal vector
% normal = [0, 0, -1];
% 
% % Unit vector for TE polarization, perpendicular to plane of incidence
% TE_direction = cross(param.beta, normal) / norm(cross(param.beta, normal));
% 
% % Check for NaN when which occurs when vector is parallel to normal
% nancheck = all(isnan(TE_direction));
% TE_direction(isnan(TE_direction)) = 0;
% 
% % When vector is parallel to normal, TE direction is defined as +y
% a_TE = TE_direction + nancheck * [0, 1, 0];
% 
% % TM polarization is defined
% a_TM = cross(a_TE, param.beta) / norm(cross(a_TE, param.beta));
% 
% % Electric field defined in surface coordinates
% p = J(1) * a_TM + J(2)* a_TE;

%% Rotation matrix

RotMat = [cos(param.theta) * cos(param.phi), -sin(param.phi), 0;...
          cos(param.theta) * sin(param.phi),  cos(param.phi), 0;...
         -sin(param.theta),                   0,              0];

p = RotMat * J;

% Magnetic fields: Hx = i*inv(mu)*(kz*Ey-ky*Ez), Hy = i*inv(mu)*(kx*Ez-kz*Ex)
h1 = 1i / param.mu_ref(iter) * (param.beta(3) * p(2) - param.beta(2) * p(3));
h2 = 1i / param.mu_ref(iter) * (param.beta(1) * p(3) - param.beta(3) * p(1));


% Electric field in Fourier space and set central harmonic to 1 which gives
% a plane wave
delta = zeros(param.num_H, 1);
delta(ceil(end / 2)) = 1;

% Incoming wave as: [Ex Ey Hx Hy]
param.s_inc = [delta * p(1); delta * p(2); delta * h1; delta * h2];
end

