function PSD(f)
f2=fftshift(fft2(f));
h=surf(log(f2.*conj(f2)/numel(f2)));
h.EdgeColor='none';
end

