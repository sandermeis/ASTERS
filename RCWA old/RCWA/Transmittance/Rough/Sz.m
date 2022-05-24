function sz = Sz(s,u)
sz=0.5*real(s{1}.*conj(1i*u{2})-s{2}.*conj(1i*u{1}));
end

