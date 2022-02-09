function input_wave = get_sinc(input_wave,device)

n_inc=sqrt(device.eps_ref*device.mu_ref);
input_wave.beta(1)=n_inc*sin(input_wave.theta)*cos(input_wave.phi);
input_wave.beta(2)=n_inc*sin(input_wave.theta)*sin(input_wave.phi);
input_wave.beta(3)=n_inc*cos(input_wave.theta);

normal=[0,0,-1];

TE_direction=cross(input_wave.beta,normal)/norm(cross(input_wave.beta,normal));
nancheck=all(isnan(TE_direction));
TE_direction(isnan(TE_direction))=0;
a_TE=TE_direction+nancheck*[0,1,0];
a_TM=cross(a_TE,input_wave.beta)/norm(cross(a_TE,input_wave.beta));

p=input_wave.pTE*a_TE+input_wave.pTM*a_TM;

p=p/norm(p);

delta = zeros(device.num_H,1);
delta(ceil(end/2))=1; %only central harmonic

h1=1i/device.mu_ref*(input_wave.beta(3)*p(2)-input_wave.beta(2)*p(3));
h2=1i/device.mu_ref*(input_wave.beta(1)*p(3)-input_wave.beta(3)*p(1));

input_wave.s_inc=[delta*p(1);delta*p(2);delta*h1;delta*h2];
end

