p=2.^(0:23);

v = VideoWriter('QR','Uncompressed AVI');
v.FrameRate=2;
open(v);

for k = 1:numel(p)
    A=readmatrix("2Doptim512K30_"+string(p(k))+".txt");
   imagesc(A)
   title("2D QR structure, " + string(p(k))+ " Iterations")
   frame = getframe(gcf);
   writeVideo(v,frame);
end

close(v);