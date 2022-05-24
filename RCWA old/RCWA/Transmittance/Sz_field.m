function sz = Sz_field(f)
[s,u]=split_xy(f);
[sx,sy]=split_xy(s);
[ux,uy]=split_xy(u);
sz=0.5*real(sx.*conj(1i*uy)-sy.*conj(1i*ux));
end

