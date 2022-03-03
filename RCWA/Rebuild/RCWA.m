function [device,layer,input_wave,Sz] = RCWA(layer,device,input_wave,wavelengthArray)
    arguments
        layer (:,1) struct
        device struct
        input_wave struct
        wavelengthArray (1,:) {mustBeNumeric,mustBeReal}
    end
    
c = onCleanup(@() progressBar(1));

if device.checkConverg
    harmArray = 1:2:device.Hmax;
else
    harmArray = device.Hmax;
end

issurf = cellfun(@(x) isa(x,"Surface"),{layer.input});
k = [layer(issurf).input];
k2 = [k.surfsize];
k3 = [k.surfres];
if ~all(k2 == k2(1))&&~all(k3 == k3(1))
    error("Mismatching layer dimensions")
end

if ~device.forceResize
    laynum = find(issurf,1);
    if laynum
        device.size_x       = layer(laynum).input.surfsize;
        device.size_y       = layer(laynum).input.surfsize;
        device.res_x        = layer(laynum).input.surfres;
        device.res_y        = layer(laynum).input.surfres;
    else
        warning("No rough layers added, forcing resize to manually entered dimensions")
    end
end

lenHarmArray = length(harmArray);

for iHarm = 1:lenHarmArray
    
    device.P            = harmArray(iHarm);
    device.Q            = harmArray(iHarm);
    device.num_H        = device.P * device.Q;
    
    device              = truncation(device);
    input_wave          = get_sinc(input_wave,device);
    layer               = build_layerstack(layer,device);
    
    lenWavelengthArray = length(wavelengthArray); 
    Sz = zeros(device.num_H,~device.calcAllRough*numel(layer)+device.calcAllRough*numel([layer.L])+2,lenWavelengthArray);

    for iWavelength = 1:lenWavelengthArray      
        progressBar(lenWavelengthArray*lenHarmArray)
        
        layer = apply_convolution(layer,device,iWavelength);
        [layer,device] = calc_K(layer,device,input_wave,iWavelength);
        Sz(:,:,iWavelength) = RCWA_transmittance(layer,device,input_wave);       
    end

R(:,iHarm) = squeeze(sum(Sz(:,end-1,:),1));

end

if device.checkConverg
    figure
    plot(R)
    xticklabels(string(wavelengthArray))
    xlabel("Wavelength (nm)")
    ylabel("Reflectance (a.u.)")
    legend(string(harmArray.^2)+' Harmonics (excl. trunc)','location','eastoutside')
end

end


function progressBar(i)

persistent iters wb tocArray;

if isempty(wb)
    wb = waitbar(0,'Initializing...');
    iters = 1;
    tic;
    tocArray = 0;
elseif iters >= i
    close(wb);
    delete(wb);
    toc
    clear wb iters tocArray
else

tocArray(end+1) = toc;

iterRemaining = i-iters;

if iters > 50
    t = tocArray(end-50:end);
else
    t = tocArray;
end

timeLeft = string(seconds(iterRemaining*mean(diff(t))),'hh:mm:ss');

waitBarDuration = iters/i;
waitBarString = {'Iteration ' + string(iters) + '/' + string(i),...
                'Estimated time remaining: ' + timeLeft};

waitbar(waitBarDuration,wb,waitBarString)

iters = iters + 1;

end

end