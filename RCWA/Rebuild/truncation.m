function device = truncation(device)
%%
M = -(device.P-1)/2:(device.P-1)/2;
N = -(device.Q-1)/2:(device.Q-1)/2;
[m,n] = meshgrid(M,N);

if device.truncateHarm
    TMAP = abs(m/((device.P-1)/2)).^(2*device.truncFactor) + abs(n/((device.Q-1)/2)).^(2*device.truncFactor);
    TMAP(isnan(TMAP))=1;
    TMAP = (TMAP <= 1);
    figure
    imagesc(TMAP);
else
    TMAP = ones(size(m));
end

% Extract Array Indices
device.tr_ind = find(TMAP(:));
device.num_H = length(device.tr_ind);

end