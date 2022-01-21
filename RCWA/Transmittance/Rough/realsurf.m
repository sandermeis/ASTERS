function gg=realsurf(eps,epsprev,Zres)
%imp=importdata("../../../ri/je882_centre_right_10um.txt");
load('RD.mat');
imp=RD{3,1}; %650 %max 487
% imp=RD{3,4}; %700 %max 1800, intermax 700
% imp=RD{3,7}; %730 %max 213
imp=-imp+max(imp(:));
N=size(imp);
d=discretize(imp,linspace(0,max(imp(:)),Zres));%25,85
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
