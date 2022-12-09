function show_fields(param, fields, direc, slice, bars)
arguments
    param
    fields
    direc {mustBeMember(direc,["x", "y", "z", "norm", "abs"])} = "norm"
    slice {mustBeMember(slice,["x", "y", "z"])} = "x"
    bars = true;
end


direc_nr = find(direc==["x", "y", "z", "norm", "abs"]);

layer = fill_layer(param);
boundaries = [layer.L];
wl = 1;
%wavelength % material % component
if direc_nr==4
    d = real(abs(fields{wl}{1})+abs(fields{wl}{2})+abs(fields{wl}{3}));
    title_text = "|E|";
elseif direc_nr==5
    d = real(abs(fields{wl}{1})+abs(fields{wl}{2})+abs(fields{wl}{3}));
    title_text = "P_{abs}";
else
    d = real(abs(fields{wl}{direc_nr}));
    title_text = "|E" + direc + "|";
end


if direc_nr~=5

plot_field_3d(d,param,layer,slice,bars,boundaries,title_text)

else
    c = 1;%2.998e8;%2*0.002654418728*256;
    for iter=1:numel(param.wavelengthArray)
        
    for i = 1:numel(layer)

        layer(i).geometry.eps = layer(i).geometry.eps_struc;
        for j = 1:numel(layer(i).permittivities)
            % Assign permittivity at specific wavelength
            eps = layer(i).permittivities{j}(param.wavelengthArray(iter)); % wl using interpolant
            layer(i).geometry.eps(layer(i).geometry.eps_struc == j) = eps;
        end
    end
    k = 1;
    for i = 1:numel(layer)
        for j = 1:numel(layer(i).L)
            for l = 1:layer(i).L(j)
                absorp(:,:,k) = -0.5 * 2*pi/param.wavelengthArray(iter) * d(:,:,k) .* imag(layer(i).geometry.eps(:,:,j));
                k = k + 1;
            end
        end
    end

    end
trapz(trapz(trapz(absorp)))
plot_field_3d(absorp,param,layer,slice,bars,boundaries,title_text)

end
% figure
% for i=1:size(newnew,3)
% imagesc(newnew(:,:,i))
% title("nr " + i)
% colorbar
% set(gca, 'clim', [cmin cmax])
% drawnow
% pause(0.3)
% end




end

function plot_field_3d(d,param,layer,slice,bars,boundaries,title_text)
            cmax = max(d, [], 'all');
            cmin = min(d, [], 'all');
    switch slice
        case {"x","y"}
            figure
            for ii = 1:size(d, 2)
                imagesc(squeeze(d(:, ii, :)))
                colorbar
                set(gca, 'clim', [cmin cmax])
                set(gca, 'YDir','normal')
                set(gca, 'XTick', [0:0.1:1] * size(d,3), 'XTickLabel', [0:0.1:1] * sum(boundaries)) % 10 ticks
                set(gca, 'YTick', [0:0.1 * param.res:param.res], 'YTickLabel', [0:0.1 * param.size / param.res:param.size / param.res] * param.res) % 20 ticks
                xlabel("z (nm)")
                ylabel(slice + " (nm)")
                title(title_text + ", " + slice + " = " + (ii-1) + " nm")
                for iii = 1:numel(layer)
                    h{iii} = text(sum([layer(1:iii).L]) - sum(layer(iii).L)/2, 0.01*param.res, layer(iii).material);
                    set(h{iii},'Rotation',90);
                end

                if bars
                    for i = 1:numel(boundaries)
                        xline(sum(boundaries(1:i)))
                    end
                end
                drawnow
                pause(0.1)
            end

        case "z"

            figure
            for ii = 1:size(d, 3)
                imagesc(squeeze(d(:, :, ii)))
                colorbar
                set(gca, 'clim', [cmin cmax])
                set(gca, 'YDir','normal')
                set(gca, 'XTick', [0:0.1 * param.res:param.res], 'XTickLabel', [0:0.1 * param.size / param.res:param.size / param.res] * param.res)
                set(gca, 'YTick', [0:0.1 * param.res:param.res], 'YTickLabel', [0:0.1 * param.size / param.res:param.size / param.res] * param.res) % 20 ticks
                xlabel("x (nm)")
                ylabel("y (nm)")
                title(title_text + ", z = " + (ii-1) + " nm")

                drawnow
                pause(0.1)
            end
    end
end