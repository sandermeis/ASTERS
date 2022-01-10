function [Stot,S] = get_power(E,mu_r,mu_r0,Kz,kz0)
Mag=abs(E{1}).^2+abs(E{2}).^2+abs(E{3}).^2;
S=real((mu_r0/mu_r)*(Kz/kz0))*Mag;
Stot=sum(S(:));
end

