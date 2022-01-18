function [Er,Etheta,Ephi]=plot_spherical(E)
Ex=E{1};
Ey=E{2};
Ez=E{3};

Er=sqrt(Ex.*conj(Ex)+Ey.*conj(Ey)+Ez.*conj(Ez));
Etheta=acos(Ez./Er); %degrees

phi=pi*(Ex<0)+pi/2*(Ex==0);
Exi=Ex;
Exi(Ex==0)=1;
Ephi=(phi+atan(Ey./Exi));
end

