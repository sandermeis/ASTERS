function layer = apply_convolution(layer, param, iter)



for i = 1:numel(layer)
    
    eps = layer(i).permittivities(param.wavelengthArray(iter));% wl
    layer(i).geometry.eps = eps * layer(i).geometry.eps_struc;
    
    if layer(i).input~=0
        if i==1
            epsprev = 1; %Air %Adapt for trn/ref
        else
            epsprev = layer(i-1).permittivities(param.wavelengthArray(iter));
        end
        layer(i).geometry.eps(layer(i).geometry.eps_struc==0) = epsprev;
        layer(i).geometry.eps = conv_mat(layer(i).geometry.eps,param.P,param.Q);
        layer(i).geometry.mu_struc = eye(param.P*param.Q);
        layer(i).geometry.mu = layer(i).geometry.mu_struc(param.tr_ind,param.tr_ind);
        layer(i).geometry.eps = layer(i).geometry.eps(param.tr_ind,param.tr_ind,:);
    end
end
end


function C = conv_mat(A,P,Q)

[Nx,Ny,zdim] = size(A);

C=zeros(P,Q,zdim);

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
                    
                    C(row,col,v) = B(p0+pfft, q0+qfft);
                end
            end
        end
    end
end

end

