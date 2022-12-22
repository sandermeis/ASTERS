function [eps_lab] = import_permittivities(materials, wavelengthArray, display, sim_nr)
arguments
    materials (1,:) cell
    wavelengthArray = []
    display = false
    sim_nr = []
end

eps_lab = cell(size(materials));

% another loop over multiple materials

for i = 1:numel(materials)
    for j = 1:numel(materials{i})
        mat = materials{i}(j);
        file_name = "src/refractive_indices/" + mat + ".csv";
        x = csvread(file_name);
        if ~isempty(wavelengthArray) && ~display
            warning('off', 'backtrace')
            % check for range of input wavelengths
            if isempty(find(x(:, 1) == wavelengthArray(1), 1))
                if wavelengthArray(1) < x(1, 1)
                    warning('Requested starting wavelength (%.2f) not in range of %s data. Extrapolating',wavelengthArray(1),mat)
                else
                    warning('Requested starting wavelength (%.2f) not in %s data. Interpolating',wavelengthArray(1),mat)
                end
            elseif isempty(find(x(:, 1) == wavelengthArray(end), 1))
                if wavelengthArray(end) > x(end, 1)
                    warning('Requested ending wavelength (%.2f) not in range of %s data. Extrapolating',wavelengthArray(end),mat)
                else
                    warning('Requested ending wavelength (%.2f) not in %s data. Interpolating',wavelengthArray(end),mat)
                end
            end
        end
        eps_lab{i}(j) = {griddedInterpolant(x(:, 1), x(:, 4))};
        ip{i}{j} = interp1(x(:,1), x(:,2:3), wavelengthArray, 'linear', 'extrap');
        %eps_lab{i} = (ip(:,1)-1i*ip(:,2)).^2;

    end
end

if display
    figure
    t = tiledlayout(1,2);
    title(t, "Material properties of layer stack, sim " + sim_nr)

    nexttile
    hold on
    k = 1;
    for i = 1:numel(materials)
        for j = 1:numel(materials{i})
            plot(wavelengthArray, ip{i}{j}(:,1), "LineWidth", 2)
            leg(k) = string(materials{i}{j});
            k = k + 1;
        end
    end
    xlabel('Wavelength (nm)')
    ylabel('Refractive index')
    xlim([wavelengthArray(1), wavelengthArray(end)])

    nexttile
    hold on
    for i = 1:numel(materials)
        for j = 1:numel(materials{i})
            plot(wavelengthArray, ip{i}{j}(:,2), "LineWidth", 2)
        end
    end
    xlabel('Wavelength (nm)')
    ylabel('Extinction coefficient')
    xlim([wavelengthArray(1), wavelengthArray(end)])
    legend(leg,'location','eastoutside');

end


end