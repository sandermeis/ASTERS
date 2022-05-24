function Q = calc_Q(eps,mu,kx,ky)
Q=1/mu*[kx*ky,eps*mu-kx^2;ky^2-eps*mu,-kx*ky];
end

