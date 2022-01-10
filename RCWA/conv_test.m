% Construct the Convolution Matrices for Materials
nh=1;

p=-nh:nh;
q=-nh:nh;

Np=length(p);
Nq=length(q);

A=ones(50);

A(10:30,10:30)=0;
subplot(2,2,1)
imagesc(abs(A));

[Nx,Ny] = size(A);
A = fftshift(fftn(A))/(Nx*Ny);  % compute fourier coefficinets of A.

subplot(2,2,2)
imagesc(abs(A));
C=zeros(Np*Nq); % convolution matrix size
p0 = 1 + floor(Nx/2);  % fftshift move the zero-order harmonic to the center of the matrix
q0 = 1 + floor(Ny/2);  % these two lines used to get to the center of A
row = 1 ; col = 1 ; % initialize row and col

for rowq = 1:Nq 
for rowp = 1:Np % scanning the row in the way of (-1, -1) (0, -1) ...
    
    
          for Q = 1:Nq % % scanning the col in the way of (-1, -1) (0, -1) ...
            for P = 1:Np
                Brow=p0+p(rowp)-p(P);
                Bcol=q0+q(rowq)-q(Q);
                sprintf('Row: %d, Col: %d',Brow,Bcol)
                C(row,col)=A(Brow,Bcol);  % get the A(p-p', q-q')
                col = col + 1;
            end
          end
        col = 1;
        row = row + 1; % move next step
end
end
subplot(2,2,3)
imagesc(abs(C));
 
 D = ifftn(A);
 subplot(2,2,4)
 imagesc(abs(D));