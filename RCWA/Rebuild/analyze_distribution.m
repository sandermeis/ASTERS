Z = rand(100,100);
M = 10;
mom = analyze(Z,M)

function mom = analyze(Z,M)
mu = mean(Z(:));
mom = zeros(1,M);
v = (Z-mu).^2;
std = sqrt(mean(v(:)));
for n = 1:M
    k = (Z - mu).^n;
    mom(n) = mean(k(:))/std^n;
end
end

