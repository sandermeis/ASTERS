function A = swap_etch_A(A)
    % 3D z-dir etch and grow
    [possible_etch, possible_grow] = pixelanalysis(A);

    where_one = find(possible_etch);
    where_zero = find(possible_grow);

    pos_one = where_one(randperm(length(where_one),1));
    pos_zero = where_zero(randperm(length(where_zero),1));

    A(pos_one)=0;
    A(pos_zero)=1;
end