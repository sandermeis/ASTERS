function [eps_lab] = import_permittivities(materials, wavelengthArray, display, sim_nr)
arguments
    materials (1,:) cell
    wavelengthArray = []
    display = false
    sim_nr = []
end

eps_lab = cell(size(materials));

for i = 1:numel(materials)
    for j = 1:numel(materials{i})
        mat = materials{i}(j);
        file_name = "src/refractive_indices/" + mat + ".csv";
        try
            x = csvread(file_name);
        catch
            error("Something went wrong when opening file " + file_name + ".")
        end

        if ~isempty(wavelengthArray) && ~display
            % check for range of input wavelengths
            if wavelengthArray(1) < x(1, 1)
                warning('Requested starting wavelength (%.2f) not in range of %s data. Extrapolating',wavelengthArray(1),mat)
            end
            if wavelengthArray(end) > x(end, 1)
                warning('Requested ending wavelength (%.2f) not in range of %s data. Extrapolating',wavelengthArray(end),mat)
            end
        end
        eps_lab{i}(j) = {griddedInterpolant(x(:, 1), x(:, 4))};
        ip{i}{j} = interp1(x(:, 1), x(:, 2:3), wavelengthArray, 'linear', 'extrap');
    end
end

% Plot permittivities
if display
    figure
    t = tiledlayout(1, 2);
    title(t, "Material properties of layer stack, sim " + sim_nr)

    nexttile
    hold on
    k = 1;
    for i = 1:numel(materials)
        for j = 1:numel(materials{i})
            plot(wavelengthArray, ip{i}{j}(:, 1), "LineWidth", 2)
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
            plot(wavelengthArray, ip{i}{j}(:, 2), "LineWidth", 2)
        end
    end
    xlabel('Wavelength (nm)')
    ylabel('Extinction coefficient')
    xlim([wavelengthArray(1), wavelengthArray(end)])
    legend(leg,'location','eastoutside');

end
end