function PQ = calc_PQ(eps,mu,Kx,Ky)

%%%                   P
%%%
%%%| Kx eps^-1 Ky     | mu - Kx eps^-1 Kx |
%%%| Ky eps^-1 Ky - mu|    - Ky eps^-1 Kx |

%%%                   Q
%%%
%%%| Kx mu^-1 Ky      | eps - Kx mu^-1 Kx |
%%%| Ky mu^-1 Ky - eps|     - Ky mu^-1 Kx |

PQ=[Kx*inv(mu)*Ky, eps-Kx*inv(mu)*Kx;
    Ky*inv(mu)*Ky-eps, -Ky*inv(mu)*Kx];

end
