%% Display fourier space
% figure
% for i=1:Nz
%     %fg=dnew2;
%     fg=A(:,:,i);%abs(fftshift(fft2(A(:,:,i))));
%     fg(round(Nx/2)+1,round(Ny/2)+1)=nan;
% imagesc(fg)
% colorbar
% drawnow
% pause(0.1)
% end

%%

%% AFM input
% rm=readmatrix("data_daan/730C_aftergrowth_centre_30um.txt");
% Nx=size(rm,1);
% Ny=size(rm,2);
% Nz=32;

% [input, Lnew, Lrecalc]=discretize_surface(rm, Nz, 0, false, false, 0, 0);
% input(input==1)=1;
% input(input==2)=0;

%%%% NEED TO PAD WITH MORE MEAT


% k_min = 20;
% k_max = 50;

% generate_QR(input)

%% Straight slab input
% Nx=256;
% Ny=256;
% Nz=32;
% 
% top=zeros(Nx,Ny,Nz);
% bot=ones(Nx,Ny,Nz);
% 
% input = cat(3,top,bot);
% 
% k_min = 35;
% k_max = 46;
% 
% iters = 4^8;
% savename = "Maarten_free_swap";
% type = "3DFree";

%% 3D input

pitch=1250;
lamb_gaas=900;
lamb0=450;

Nx=64;
Ny=64;
Nz=32;

input = zeros(Nx,Ny,Nz);

input(randperm(Nx*Ny*Nz,round(Nx*Ny*Nz/2)))=1;

ngaas=3.5975;%1.42;

k_min = (pitch/lamb0);
k_max = (pitch/lamb_gaas)*ngaas;

num_iter = 17;
iv = (5:num_iter);


p2size=round(log2(Nx*Ny*Nz));

sizes=2.^p2size./(2.^iv);

intervals = 2.^iv;

samp = cumsum(intervals);

iters = sum(intervals);

savename = "Maarten_free_swap_med3";
type = "3DFree";

%% 2D input
% pitch=1350;
% lamb_gaas=900;
% lamb0=450;
% 
% Nx=51;
% Ny=51;
% 
% input = zeros(Nx,Ny);
% 
% input(randperm(Nx*Ny,round(Nx*Ny/2)))=1;
% 
% k0=2*pi/lamb0;
% 
% ngaas=3.5975;
% kgaas=ngaas*2*pi/lamb_gaas;
% 
% spacing = 2*pi/pitch;
% 
% k_min = k0/spacing;
% k_max = kgaas/spacing;
% 
% iters = 2^23;
% savename = "2Doptim51K3";
% type = "2D";
%%
generate_QR(input, samp, intervals, sizes, k_min, k_max, type, savename)

%%
function generate_QR(A, samp, intervals, sizes, k_min, k_max, type, savename)
arguments
    A {isnumeric}
    samp {isnumeric}
    intervals {isnumeric}
    sizes {isnumeric}
    k_min (1,1) {isnumeric}
    k_max (1,1) {isnumeric}
    type {isStringScalar}
    savename {isStringScalar}
end

tic

iters = samp(end);

[Nx,Ny,Nz]=size(A);

B=zeros(Nx,Ny,Nz);

[kx,ky,kz]=ndgrid(-(Nx-1)/2:(Nx-1)/2,-(Ny-1)/2:(Ny-1)/2,-(Nz-1)/2:(Nz-1)/2);
k2=sqrt(kx.^2+ky.^2);

B(k2>k_min & k2<k_max)=1;
%%
% [sx,sy]=size(A);
% 
% B=zeros(sx,sy);
% %%
% [kx,ky]=ndgrid(-(sx-1)/2:(sx-1)/2,-(sy-1)/2:(sy-1)/2);
% k2=sqrt(kx.^2+ky.^2);


%B(k2>k_min&k2<k_max)=1;
num_target=sum(B,'all');


%% Optimization constants
y=0.5;
M=y*(1-y)*Nx^4/num_target;

%zeros(6*H+1);
%A(c2(1:round(N2*y)))=1;
%A((2*H+1):(4*H+1),(2*H+1):(4*H+1))=A;

%zeros(6*H+1);
%B((2*H+1):(4*H+1),(2*H+1):(4*H+1))=B;
%%
    f_diff = (abs(fftshift(fftn(A)))-M).^2;

    % Transformation
    %f_diff = abs(fftshift(fftn(A_try)));

    % Compare with target
    T = sum(f_diff(logical(B)),'all');


successes = 0;
fails = 0;
consecutive_fails = 0;

% Display progress
tiledlayout("flow")
nexttile
imagesc(sum(A,3))
title("Swap 0")
drawnow

%% Gaussian ordering

% Attenuating categories
% cats=12;
% s=Nx/2*1./kron(1:cats,ones(1,iters/cats));
% h=Nz/8*1./kron(1:4,ones(1,iters/4));


% Straight

% x-y size
%s=Nx/16*ones(1,iters);

% height
%h=Nz/8*ones(1,iters);
%%

writematrix(A,"generated_surfaces/" + savename + ".txt")


% Write intervals
%samp = 2.^(1:num_iter);
p=1;



for i = 1:iters

    switch type
        case "2D"
            A_try = shift_A(A);
        case "3DBlock"
            %for i2 = 1:200
                A_try = swap_etch_A(A);
            %end
        case "3DFree"
            %for i2 = 1:200
                A_try = swap_random_A(A,sizes(p));
            %end
        case "3DGauss"
            A_try = dent_A(A,s(i),h(i));
    end
    %A_try((2*H+1):(4*H+1),(2*H+1):(4*H+1)) = swap_A(A((2*H+1):(4*H+1),(2*H+1):(4*H+1)));
    
    %% Optimization
    f_diff = (abs(fftshift(fftn(A_try)))-M).^2;

    % Transformation
    %f_diff = abs(fftshift(fftn(A_try)));

    % Compare with target
    T_try = sum(f_diff(logical(B)),'all');

    if T_try < T
        A = A_try;
        T = T_try;
        consecutive_fails = 0;
        successes = successes + 1;
    else
        fails = fails + 1;
        consecutive_fails = consecutive_fails + 1;
    end
    if i == samp(p)
        nexttile
        imagesc(sum(A,3))
        title("amt: "+intervals(p)+"sz: "+sizes(p)+"itr: "+i)
        drawnow
        disp("Iterations: " + i)
        disp("-----------------------")
        disp("Successes: " + successes)
        disp("Fails: " + fails)
        disp("Consecutive fails: " + consecutive_fails)
        disp("-----------------------")
        writematrix(A, "generated_surfaces/" + savename + "_" + i + ".txt")
        p = p + 1;
    end
    if consecutive_fails == 200000
        disp("Converged")
        break
    end
end

toc

end
% disp("shifting")
% iters = 50000;
% fails=0;
% successes=0;
% for i = 1:iters
%     A_try = shift_A(A);
%     f_diff = (abs(fftshift(fft2(A_try)))-M).^2;
%     T_try = sum(f_diff(logical(B)),'all');
% 
%     if T_try<T
%         A = A_try;
%         T = T_try;
%         disp("success")
%         fails=0;
%         successes=successes+1;
%     else
%         fails=fails+1;
%     end
%     if mod(successes,100)==0
%         imagesc(A)
%         title("shift "+i)
%         drawnow
%     end
%     if fails==10000
%         disp("converged")
%         break
%     end
%     disp(i)
% end
%%

% bbb=fftshift(fftn(A));
% bbb(round((end+1)/2),round((end+1)/2))=NaN;
% 
% figure
% tiledlayout("flow")
% 
% nexttile
% imagesc(abs(bbb))
% 
% nexttile
% surf(1:size(bbb,1),1:size(bbb,2),real(bbb),imag(bbb),'EdgeColor','none','FaceColor','interp')
% colorbar
% 
% nexttile
% imagesc(abs(bbb))
% colorbar

% figure
% h=slice(1:Nx,1:Ny,1:(Nz+Nz/8),A,[],[],1:5:Nz);
% set(h,'edgecolor','none')

% xlabel('x')
% ylabel('y')
% zlabel('z')

