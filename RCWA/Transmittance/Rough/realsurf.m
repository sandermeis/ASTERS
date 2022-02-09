function gg=realsurf(eps,epsprev,Zres)
%imp=importdata("../../../ri/je882_centre_right_10um.txt");
load('RD.mat','RD');
reverse=false;
sf=730;

switch sf
    case 650
        imp=RD{3,1}; %650 %max 487
    case 700
        imp=RD{3,4}; %700 %max 1800, intermax 700
    case 730
        imp=RD{3,7}; %730 %max 213
end
imp=imp-min(imp(:));

N=size(imp);
d=discretize(imp,linspace(0,max(imp(:)),Zres));%25,85
% figure
% h2=surf(d);
% h2.EdgeColor='none';
%
gg=zeros(N(1),N(2),Zres);

if reverse
    for i=1:Zres
        gg(:,:,Zres-i+1)=(d>=i);
    end
else
    for i=1:Zres
        gg(:,:,i)=(d>=i);
    end
    
end

gg=gg*eps;
gg(gg==0)=epsprev;
end
