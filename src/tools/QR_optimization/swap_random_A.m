function A = swap_random_A(A, amount)
    % 3D random block etch and grow
    %[possible_etch, possible_grow] = pixelanalysis(A);

    where_one = find(A);
    where_zero = find(~A);

    pos_one = where_one(randperm(length(where_one),amount));
    pos_zero = where_zero(randperm(length(where_zero),amount));

    A(pos_one)=0;
    A(pos_zero)=1;
end