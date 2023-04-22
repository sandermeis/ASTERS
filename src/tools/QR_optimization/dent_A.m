function A_out = dent_A(A,s,h)
% Adding and/or removing gaussian surfaces
%%
A_out=A;
[Nx,Ny,Nz]=size(A_out);
sx=s;
sy=s;

A_flat=sum(A_out,3);

%% retry if fail
for tries=1:100
[pos_growx,pos_growy] = ind2sub([Nx,Ny],randperm(Nx*Ny,1));
[pos_shrinkx,pos_shrinky] = ind2sub([Nx,Ny],randperm(Nx*Ny,1));

gpos_grow = gauss2d(Nx, pos_growx, pos_growy, sx, sy, h);
gpos_shrink = -gauss2d(Nx, pos_shrinkx, pos_shrinky, sx, sy, h);


    %disp("Try: " + tries)
    % pbc
if all((A_flat+gpos_grow)<=Nz,'all')&&all((A_flat+gpos_shrink)>0,'all')
    %do it
    %disp("Success")
    [x,y]=ndgrid(1:Nx,1:Ny);
    for i=1:Nx*Ny

    A_out(x(i),y(i),(A_flat(x(i),y(i))+1):(A_flat(x(i),y(i))+gpos_grow(x(i),y(i))))=1;
    A_out(x(i),y(i),((A_flat(x(i),y(i))+1+gpos_shrink(x(i),y(i))):A_flat(x(i),y(i))))=0;
    end
% tiledlayout("flow")
% nexttile
% contour(A_flat)
% nexttile
% contour(sum(A_out,3))
    break
end
end
end