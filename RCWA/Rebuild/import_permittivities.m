function [eps_lab] = import_permittivities(materials)%, wl)
arguments
    materials (1,:) cell
%     wl (1,:) cell
end

eps_lab = cell(size(materials));

for i=1:numel(materials)
%     wavelengthArray = wl{i};
    mat = materials{i};
    file_name = "refractive_indices/"+mat+".csv";
    x = csvread(file_name);

%     warning('off','backtrace')
%     % check for range of input wavelengths
%     if isempty(find(x(:,1)==wavelengthArray(1),1))
%         if wavelengthArray(1)<x(:,1)
%             warning('Requested starting wavelength (%.2f) not in range of %s data. Extrapolating',wavelengthArray(1),mat)
%         else
%             warning('Requested starting wavelength (%.2f) not in %s data. Interpolating',wavelengthArray(1),mat)
%         end
%     elseif isempty(find(x(:,1)==wavelengthArray(end),1))
%         if wavelengthArray(end)>x(:,1)
%             warning('Requested ending wavelength (%.2f) not in range of %s data. Extrapolating',wavelengthArray(end),mat)
%         else
%             warning('Requested ending wavelength (%.2f) not in %s data. Interpolating',wavelengthArray(1),mat)
%         end
%     end
    eps_lab{i} = griddedInterpolant(x(:,1),x(:,4));
    %ip = interp1(x(:,1), x(:,2:3), wavelengthArray, 'linear', 'extrap');
    %eps_lab{i} = (ip(:,1)-1i*ip(:,2)).^2;

end

end

