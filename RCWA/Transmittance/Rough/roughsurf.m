function [gg]=roughsurf(eps,epsprev,Xres,Yres,Zres,rr)
%%
rng(1337)
N = [Xres Yres];
F = rr;
[X,Y] = ndgrid(1:N(1),1:N(2));
i = min(X-1,N(1)-(X-1));
j = min(Y-1,N(2)-(Y-1));

% qr=
% qL=
% qs=
% H=
H = exp(-.5*(i.^2+j.^2)/F^2);
Z = real(ifft2(H.*fft2(randn(N))));

%normalize
data=Z/max(max(abs(Z)));
% figure
% h1=surf(data);
% h1.EdgeColor='none';
%%
% rmsr=sqrt(sum(sum(data.^2))/numel(data))
% PSD(data)
d=discretize(data,linspace(-1,1,Zres));
% figure
% h2=surf(d);
% h2.EdgeColor='none';
%%
gg=zeros(N(1),N(2),Zres);
for i=1:Zres
gg(:,:,i)=i<=d;
end
gg=gg*eps;
gg(gg==0)=epsprev;
end