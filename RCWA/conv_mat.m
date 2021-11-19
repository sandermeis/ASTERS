function C=conv_mat(A,P,Q)

[Nx,Ny] = size(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Number of Harmonics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nharmonics = P * Q;

p = [-floor(P/2):+floor(P/2)];
q = [-floor(Q/2):+floor(Q/2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Fourier Coefficients of A
%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = fftshift(fftn(A)) / (Nx*Ny);

%%%%%%%%%%%%%%%%
% Array Indicies
%%%%%%%%%%%%%%%%
p0 = floor(Nx/2) + 1; %Added 1
q0 = floor(Ny/2) + 1; %because of

%%%%%%%%%%%%%%%%%%%%%
% Loop Over Harmonics
%%%%%%%%%%%%%%%%%%%%%

for qrow = 1:Q
    for prow = 1:P
        row = (qrow-1)*P + prow;
        for qcol = 1:Q
            for pcol = 1:P
                col =  (qcol - 1)*P + pcol;
                pfft = p(prow) - p(pcol);
                qfft = q(qrow) - q(qcol);
                
                C(row,col) = B(p0+pfft, q0+qfft);
            end %pcol
        end %qcol
    end %prow
end %qrow

end