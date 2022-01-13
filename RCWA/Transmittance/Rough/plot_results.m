function [h] = plot_results(lam0_r,n,J,Rtot,Ttot)

h = figure;
lab1 = lam0_r(1);
lab2 = lam0_r(end);

N=size(J,2);

x_grid=repmat(lam0_r',1,N+2);

J(:,N+1)=Rtot;
J(:,N+2)=Ttot;

n(N+1)="R";
n(N+2)="T";

area(x_grid,J)
% hold on
% plot(lam0_r,real(eps_lab{1})/max(real(eps_lab{1})))
% plot(lam0_r,imag(eps_lab{1})/max(real(eps_lab{1})))
legend(n)
%ylim([-0.5,1.5])
xlim([lab1-0.1*(lab2-lab1),lab2+0.1*(lab2-lab1)])
end

