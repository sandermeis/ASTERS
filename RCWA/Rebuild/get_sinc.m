function param = get_sinc(param, iter)

n_inc = sqrt(param.eps_ref(iter)*param.mu_ref(iter));
param.beta(1) = n_inc * sin(param.theta) * cos(param.phi);
param.beta(2) = n_inc * sin(param.theta) * sin(param.phi);
param.beta(3) = n_inc * cos(param.theta);

normal = [0, 0, -1];

TE_direction = cross(param.beta, normal) / norm(cross(param.beta, normal));
nancheck = all(isnan(TE_direction));
TE_direction(isnan(TE_direction)) = 0;
a_TE = TE_direction + nancheck * [0,1,0];
a_TM = cross(a_TE, param.beta) / norm(cross(a_TE, param.beta));

p = param.pTE * a_TE + param.pTM * a_TM;

p = p / norm(p);

delta = zeros(param.num_H,1);
delta(ceil(end/2)) = 1; %only central harmonic

h1 = 1i / param.mu_ref(iter) * (param.beta(3) * p(2) - param.beta(2) * p(3));
h2 = 1i / param.mu_ref(iter) * (param.beta(1) * p(3) - param.beta(3) * p(1));

param.s_inc = [delta * p(1); delta * p(2); delta * h1; delta * h2];
end

