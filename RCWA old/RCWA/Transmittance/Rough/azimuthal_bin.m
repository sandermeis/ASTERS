function thbin=azimuthal_bin(TH,R,th_range)

for i=1:length(th_range)-1
    thbin(i)=mean(R(TH>=th_range(i) & TH<th_range(i+1)));  
end
end

