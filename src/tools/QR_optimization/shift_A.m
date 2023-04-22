function A = shift_A(A)
    % 2D binary shift blocks
    pos = randperm(numel(A),1);
    [sy,sx]=size(A);

    positions={[0,-1],[0,1],[-1,0],[1,0]};
    order=randperm(4,4);
    eligible=[];
    for i=1:4
        test_pos(i) = toPBC(sx,sy,pos,positions{order(i)});
        if A(pos)~=A(test_pos(i))
            eligible(i)=1;
        end
    end
    if ~isempty(eligible)
    f_e = find(eligible);
    tp=test_pos(f_e(randperm(numel(f_e),1)));

    temp=A(pos);
    A(pos)=A(tp);
    A(tp)=temp;
    end

end

function lind_pbc = toPBC(sx,sy,pos,shift)
        r = mod(pos-1,sx)+1+shift(1);
        c = ceil(pos./sx)+shift(2);
        r_pbc = mod(r-1,sx)+1;
        c_pbc = mod(c-1,sy)+1;
        lind_pbc = (c_pbc-1).*sx+r_pbc;
end

