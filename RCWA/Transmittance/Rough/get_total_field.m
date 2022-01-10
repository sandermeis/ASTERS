function [s_f,s_b,u_f,u_b] = get_total_field(Kx,Ky,eps_r,mu_r,W,V,lab,cp,cm,z)
[s_f{1},s_f{2}]=split_xy(W*expm(-lab*z)*cp);
[s_b{1},s_b{2}]=split_xy(W*expm(lab*z)*cm);
[u_f{1},u_f{2}]=split_xy(-V*expm(-lab*z)*cp);
[u_b{1},u_b{2}]=split_xy(V*expm(lab*z)*cm);

s_f{3}=-1i*inv(eps_r)*(Kx*u_f{2}-Ky*u_f{1});
s_b{3}=-1i*inv(eps_r)*(Kx*u_b{2}-Ky*u_b{1});
u_f{3}=-1i*inv(mu_r)*(Kx*s_f{2}-Ky*s_f{1});
u_b{3}=-1i*inv(mu_r)*(Kx*s_b{2}-Ky*s_b{1});
end

