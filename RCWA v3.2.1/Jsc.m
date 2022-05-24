function jsc = Jsc(eqe, wavelengthArray)

    hcq = 1986.4458571489287/1.602176634; %hc/q: J*nm/C
    ssi = readmatrix("Sun_spectrum_AM15g.csv");
    ss = interp1(ssi(:,1),ssi(:,2),wavelengthArray,'linear','extrap');
for i=1:size(eqe,1)
    jsc(i) = trapz(wavelengthArray,eqe(i,:).*ss.*wavelengthArray)/hcq*0.1; %C * Int[frac_ph * ph/m^2/nm]nm = I/m^2 = 10mA/cm^2
end

end