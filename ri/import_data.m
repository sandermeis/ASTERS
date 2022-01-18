im=importdata("JE882_reflectance.qe");
figure
leg=[];
for i=1:length(im.colheaders)
    if im.colheaders{i}=="Wavelength [nm]"
        plot(im.data(:,i),im.data(:,i+1));
        if isempty(leg)
        leg=[string(im.colheaders{i+1})];
        else
            leg=[leg,string(im.colheaders{i+1})];
        end
        hold on
    end
end
legend(leg)

im2=importdata("je882_centre_right_10um.txt");
figure
imagesc(im2*1e9)
colorbar()