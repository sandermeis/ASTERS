function eps_lab=get_lab(lab1,lab2)
%% import and check data

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS"];

for i0=1:length(materials)
    FileName=strcat(["../ri/refractive indices.xlsx - "],materials{i0},".csv");
    x{i0} = csvread(FileName(1,1)); %x{i0} is tab of material i0
    
    % check for range of input wavelengths
    if isempty(find(x{i0}(:,1)==lab1,1))
        error('Error. \nRequested starting wavelength (%s) not in range of %s data.',lab1,materials{i0})
    elseif isempty(find(x{i0}(:,1)==lab2,1))
        error('Error. \nRequested ending wavelength (%s) not in range of %s data.',lab2,materials{i0})
    end
    
    ind1=find(x{1,i0}(:,1)==lab1);
    ind2=find(x{1,i0}(:,1)==lab2);
    
    eps_lab{i0}=(x{1,i0}(ind1:ind2,2)-1i*x{1,i0}(ind1:ind2,3)).^2; %(second row-i*third row)^2

end

end

