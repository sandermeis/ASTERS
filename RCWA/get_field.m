function E = get_field(Kx,Ky,Kz,W,lab,cp,cm,z)
[E{1},E{2}]=split_xy(W*expm(-lab*z)*cp+W*expm(lab*z)*cm);
E{3}=-inv(Kz)*(Kx*E{1}+Ky*E{2});
end

