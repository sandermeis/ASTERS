function eps_lab = get_lab(lam0_r)

lab1=lam0_r(1);
lab2=lam0_r(end);
% dlab=lam0_r(2)-lam0_r(1);
% import and check data

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS","GaP","AlInP"];

for i=1:length(materials)
    FileName=strcat(["../../../ri/refractive indices.xlsx - "],materials{i},".csv");
    x{i} = csvread(FileName(1,1)); %x{i0} is tab of material i0
    warning('off','backtrace')
    % check for range of input wavelengths
    if isempty(find(x{i}(:,1)==lab1,1))
        if lab1<x{i}(:,1)
            warning('Requested starting wavelength (%.2f) not in range of %s data. Extrapolating',lab1,materials{i})
        else
            warning('Requested starting wavelength (%.2f) not in %s data. Interpolating',lab1,materials{i})
        end
    elseif isempty(find(x{i}(:,1)==lab2,1))
        if lab2>x{i}(:,1)
            warning('Requested ending wavelength (%.2f) not in range of %s data. Extrapolating',lab2,materials{i})
        else
            warning('Requested ending wavelength (%.2f) not in %s data. Interpolating',lab1,materials{i})
        end
    end
    
%     ind1=find(x{1,i}(:,1)==lab1);
%     ind2=find(x{1,i}(:,1)==lab2);
    ip=interp1(x{i}(:,1),x{i}(:,2:3),lam0_r,'linear','extrap');
    eps_lab{i}=(ip(:,1)-1i*ip(:,2)).^2;
%     eps_lab2{i}=(x{1,i}(ind1:dlab:ind2,2)-1i*x{1,i}(ind1:dlab:ind2,3)).^2; %(second row-i*third row)^2
     
end

end

