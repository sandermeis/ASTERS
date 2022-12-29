function jsc_out = jsc(eqe, wavelengthArray)

hc = 6.62607015e-34 * 299792458; %hc/q: V*nm
q = 1.602176634e-19;
ssi = readmatrix("src/data/Sun_spectrum_AM15g.csv");
if numel(wavelengthArray) == 1
    jsc_out = 0;
    warning("Jsc undefined for 1 wavelength input")
else
    for i = 1:size(eqe,1)
        temp_eqe = eqe(i, :);
        where_nan = isnan(temp_eqe);
        if any(where_nan)
            warning("NaN(s) removed from EQE")
            temp_eqe = temp_eqe(~where_nan);
        end
        temp_wavelengthArray = wavelengthArray(~where_nan);

        ss = interp1(ssi(:, 1), ssi(:, 2), temp_wavelengthArray, 'linear', 'extrap');

        % C * Int[frac_ph * ph/m^2/nm]nm = I/m^2 = 10mA/cm^2
        jsc_out(i) = 0.1 * q * trapz(temp_wavelengthArray, temp_eqe .* ss .* temp_wavelengthArray) * 1e-9 / hc;
    end
end

end