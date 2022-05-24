function [a,b] = split_xy(arr)
n=length(arr);
a=arr(1:n/2);
b=arr(n/2+1:end);
end

