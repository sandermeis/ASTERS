function Cn=conv_mat(A,P,Q)

[Nx,Ny,zdim] = size(A);

for v=1:zdim

p = -floor(P/2):floor(P/2);
q = -floor(Q/2):floor(Q/2);

B = fftshift(fft2(A))/(Nx*Ny);

p0 = floor(Nx/2) + 1;
q0 = floor(Ny/2) + 1;

for qrow = 1:Q
    for prow = 1:P
        row = (qrow-1)*P + prow;
        for qcol = 1:Q
            for pcol = 1:P
                col =  (qcol - 1)*P + pcol;
                pfft = p(prow) - p(pcol);
                qfft = q(qrow) - q(qcol);
                
                C(row,col) = B(p0+pfft, q0+qfft);
            end
        end
    end
end

Cn(:,:,v)=C;

end

end