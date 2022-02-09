function layer = import_permittivities(layer,lam0_r,plot_permfig)
    arguments
        layer (:,1) struct
        lam0_r (1,:) {mustBeNumeric,mustBeReal}
        plot_permfig logical  = true
    end
    
lab1=lam0_r(1);
lab2=lam0_r(end);

for i=1:numel(layer)
    mat = layer(i).material;
    file_name = "../../ri/refractive indices.xlsx - "+mat+".csv";
    x{i} = csvread(file_name(1,1));
    
    warning('off','backtrace')
    % check for range of input wavelengths
    if isempty(find(x{i}(:,1)==lab1,1))
        if lab1<x{i}(:,1)
            warning('Requested starting wavelength (%.2f) not in range of %s data. Extrapolating',lab1,mat)
        else
            warning('Requested starting wavelength (%.2f) not in %s data. Interpolating',lab1,mat)
        end
    elseif isempty(find(x{i}(:,1)==lab2,1))
        if lab2>x{i}(:,1)
            warning('Requested ending wavelength (%.2f) not in range of %s data. Extrapolating',lab2,mat)
        else
            warning('Requested ending wavelength (%.2f) not in %s data. Interpolating',lab1,mat)
        end
    end
    

    ip{i}=interp1(x{i}(:,1),x{i}(:,2:3),lam0_r,'linear','extrap');
    eps_lab{i}=(ip{i}(:,1)-1i*ip{i}(:,2)).^2;
  
end

[layer.permittivities]=deal(eps_lab{:});

if plot_permfig
    t = tiledlayout(1,2);
    title(t,'Material properties of layer stack')
    nexttile
    hold on
    for i=1:length(ip)
        plot(lam0_r,ip{i}(:,1))
    end
    xlabel('Wavelength (nm)')
    ylabel('Refractive index')
    xlim([lab1,lab2])
    
    nexttile
    hold on
    for i=1:length(eps_lab)
        plot(lam0_r,ip{i}(:,2))
    end
    xlabel('Wavelength (nm)')
    ylabel('Extinction coefficient')
    xlim([lab1,lab2])
    legend(string({layer.material}),'location','eastoutside');
end

end

