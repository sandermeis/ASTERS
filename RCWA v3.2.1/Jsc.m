function jsc = Jsc(eqe, wavelengthArray)

% for nan
% for 1 wavelength
hcq = 1986.4458571489287/1.602176634; %hc/q: V*nm
ssi = readmatrix("Sun_spectrum_AM15g.csv");
for i=1:size(eqe,1)
    if numel(wavelengthArray)==1
        jsc(i)=0;
        if i==1
            warning("Jsc undefined for 1 wavelength input")
        end
    else

        temp_eqe=eqe(i,:);
        where_nan=isnan(temp_eqe);
        if any(where_nan)
            warning("NaN(s) removed from EQE")
            temp_eqe=temp_eqe(~where_nan);
        end
        temp_wavelengthArray=wavelengthArray(~where_nan);

        ss = interp1(ssi(:,1),ssi(:,2),temp_wavelengthArray,'linear','extrap');


        jsc(i) = trapz(temp_wavelengthArray,temp_eqe.*ss.*temp_wavelengthArray)/hcq*0.1; %C * Int[frac_ph * ph/m^2/nm]nm = I/m^2 = 10mA/cm^2
    end
end

end