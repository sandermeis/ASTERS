function lind_pbc = toPBC(NX,NXnew,NYnew,ind)
r = mod(ind-1,NX)+1;
c = ceil(ind./NX);
r_pbc = mod(r-1,NXnew)+1;
c_pbc = mod(c-1,NYnew)+1;
lind_pbc = (c_pbc-1).*NXnew+r_pbc;
end
