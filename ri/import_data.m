im_650=importdata("650C_totaldiffusereflectance_inverted_noARC_QE.txt");
im_650_ARC=importdata("650C_totaldiffusereflectance_inverted_ARC_centre_QE.txt");
im_650_AG=importdata("650C_as-grown.xlsx");

im_700=importdata("700C_totaldiffusereflectance_inverted_noARC_QE.txt");
im_700_ARC=importdata("700C_totaldiffusereflectance_inverted_ARC_QE.txt");
im_700_AG=importdata("700C_as-grown.xlsx");

im_730=importdata("730C_totaldiffusereflectance_inverted_noARC_QE.txt");
im_730_ARC=importdata("730C_totaldiffusereflectance_inverted_ARC_QE.txt");
im_730_AG=importdata("730C_as-grown.xlsx");

im_650_STM=1e9*importdata("650C_aftergrowth_centre_30um.txt");
im_700_STM=1e9*importdata("700C_aftergrowth_centre_30um.txt");
im_730_STM=1e9*importdata("730C_aftergrowth_centre_30um.txt");

im_700_STM_fixed=im_700_STM;
im_700_STM_fixed(1:118,:)=im_700_STM_fixed(1:118,:)-165;
im_700_STM_fixed(119,1:232)=im_700_STM_fixed(119,1:232)-165;

%still some funky stuff going on
%plotfig(im_650,im_700,im_730)

% RD{1,1}=im_650.textdata(2,3);
% RD{1,2}=im_650.textdata(2,4);
% RD{1,3}=im_650.textdata(2,6);
% RD{1,4}=im_700.textdata(2,1);
% RD{1,5}=im_700.textdata(2,2);
% RD{1,6}=im_700.textdata(2,4);
% RD{1,7}=im_730.textdata(2,3);
% RD{1,8}=im_730.textdata(2,4);
% RD{1,9}=im_730.textdata(2,6);
% 
% RD{2,1}=im_650.data(1:end,3);
% RD{2,2}=im_650.data(1:end,4);
% RD{2,3}=im_650.data(1:end,6);
% RD{2,4}=im_700.data(1:end,1);
% RD{2,5}=im_700.data(1:end,2);
% RD{2,6}=im_700.data(1:end,4);
% RD{2,7}=im_730.data(1:end,3);
% RD{2,8}=im_730.data(1:end,4);
% RD{2,9}=im_730.data(1:end,6);
% 
% RD{3,1}=im_650_STM;
% RD{3,4}=im_700_STM_fixed;
% RD{3,7}=im_730_STM;
% 
% save('RD.mat','RD')
% fsc=importdata("TMM_res.xlsx");
% area(fsc.data.FullSolarCell(1:701,1),fsc.data.FullSolarCell(1:701,2:10));
% legend(fsc.textdata.FullSolarCell(1,2:10))

fscf=importdata("TMM_res_fix.xlsx");
figure;
area(fscf.data(1:701,1),fscf.data(1:701,2:10));
legend(fscf.textdata(2,2:10))

imagesc(100*([J,Rtot',Ttot']-fscf.data(1:701,2:10)));
c = colorbar;  
c.Ruler.TickLabelFormat='%g%%';
title("Difference in TMM and RCWA Absorption")
xticklabels([n,"R","T"])
ylabel("Wavelength (Nm)")

function plotfig(varargin)
switch nargin
    case 1
        im=varargin{1};
        if ~isstruct(varargin{1})
            
            figure
            imagesc(im*1e9)
            colorbar()
        else
            figure
            leg=[];
            for i=1:(length(im.colheaders))
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
        end
    otherwise
        figure
        leg=[];
        hold on
        for j=1:nargin
            im=varargin{j};
            kz=3-2*(j==2);
            for i=kz:kz%(length(im.colheaders))
                if im.colheaders{i}=="Wavelength [nm]"
                    plot(im.data(:,i),im.data(:,i+3)./im.data(:,i+1),'Marker', 'o');
                    if isempty(leg)
                        leg=[string(im.colheaders{i+1})];
                    else
                        leg=[leg,string(im.colheaders{i+1})];
                    end
                    
                end
            end
            
        end
        legend(leg)
end
end
