function [E_inc1,E_inc2] = get_Einc(phi,theta,pte,ptm,num_H)

    n_inc=1; % needs generalization
    k_inc = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]';%*k_0*n_inc;
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
    
    E_inc1 = delta*p(1);
    E_inc2 = delta*p(2);
end

