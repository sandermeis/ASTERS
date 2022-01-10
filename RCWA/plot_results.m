function plot_results(lam0_r,n,J,Rtot,Ttot)

lab1 = lam0_r(1);
lab2 = lam0_r(end);

for i=1:length(J)
Abs(:,i)=J{i};
pll(:,i)=lam0_r';
leg(i)=n(i);
if i==length(J)
Abs(:,i+1)=Rtot;
pll(:,i+1)=lam0_r';
leg(i+1)="R";

Abs(:,i+2)=Ttot;
pll(:,i+2)=lam0_r';
leg(i+2)="T";
end
end
area(pll,Abs)
% hold on
% plot(lam0_r,real(eps_lab{1})/max(real(eps_lab{1})))
% plot(lam0_r,imag(eps_lab{1})/max(real(eps_lab{1})))
legend(leg)
%ylim([-0.5,1.5])
xlim([lab1-0.1*(lab2-lab1),lab2+0.1*(lab2-lab1)])
end

