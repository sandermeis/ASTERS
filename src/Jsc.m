function jsc = Jsc(eqe, wavelengthArray)

hcq = 1986.4458571489287/1.602176634; %hc/q: V*nm
ssi = readmatrix("Sun_spectrum_AM15g.csv");
if numel(wavelengthArray) == 1
    jsc = 0;
    warning("Jsc undefined for 1 wavelength input")
else
    for i = 1:size(eqe,1)
        temp_eqe = eqe(i,:);
        where_nan = isnan(temp_eqe);
        if any(where_nan)
            warning("NaN(s) removed from EQE")
            temp_eqe = temp_eqe(~where_nan);
        end
        temp_wavelengthArray = wavelengthArray(~where_nan);

        ss = interp1(ssi(:, 1), ssi(:, 2), temp_wavelengthArray, 'linear', 'extrap');

        % C * Int[frac_ph * ph/m^2/nm]nm = I/m^2 = 10mA/cm^2
        jsc(i) = trapz(temp_wavelengthArray, temp_eqe .* ss .* temp_wavelengthArray) / hcq * 0.1;
    end
end

end