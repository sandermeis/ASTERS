function param = truncation(param)
%%
for i=1:numel(param)
    M = -(param(i).P-1)/2:(param(i).P-1)/2;
    N = -(param(i).Q-1)/2:(param(i).Q-1)/2;
    [m,n] = meshgrid(M,N);

    if param(i).truncateHarm
        TMAP = abs(m/((param(i).P-1)/2)).^(2*param(i).truncFactor) + abs(n/((param(i).Q-1)/2)).^(2*param(i).truncFactor);
        TMAP(isnan(TMAP))=1;
        TMAP = (TMAP <= 1);
        if param(i).truncfig
            figure
            imagesc(TMAP);
        end
    else
        TMAP = ones(size(m));
    end

    % Extract Array Indices
    param(i).tr_ind = find(TMAP(:));
    param(i).num_H = length(param(i).tr_ind);
end


end