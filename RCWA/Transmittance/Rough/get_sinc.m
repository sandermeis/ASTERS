function [s_inc] = get_sinc(k_0,eps,mu,phi,theta,pte,ptm,num_H)

    k_inc = sqrt(eps*mu)*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]';
    normal=[0;0;-1];
    
    TE_direction=cross(k_inc,normal)/norm(cross(k_inc,normal));
    nancheck=all(isnan(TE_direction));
    TE_direction(isnan(TE_direction))=0;
    a_TE=TE_direction+nancheck*[0,1,0]';
    a_TM=cross(a_TE,k_inc)/norm(cross(a_TE,k_inc));
    
    p=pte*a_TE+ptm*a_TM;
    
    p=p/norm(p);
    
    delta=zeros(num_H,1);
    delta(ceil(end/2))=1; %only central harmonic
    
    h1=1i/mu*(k_inc(3)*p(2)-k_inc(2)*p(3));
    h2=1i/mu*(k_inc(1)*p(3)-k_inc(3)*p(1));
    
    s_inc=[delta*p(1);delta*p(2);delta*h1;delta*h2];
end

