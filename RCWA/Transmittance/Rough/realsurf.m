function gg=realsurf(eps,epsprev,Zres)
imp=importdata("../../../ri/je882_centre_right_10um.txt");
imp=imp*1e9; %to nm
Zres=25;
N=size(imp);
d=discretize(imp,linspace(0,425,Zres));%25,85
% % figure
% % h2=surf(d);
% % h2.EdgeColor='none';
gg=zeros(N(1),N(2),Zres);
for i=1:Zres
gg(:,:,i)=i<=d;
end
gg=gg*eps;
gg(gg==0)=epsprev;
end
