function [device,layer,input_wave,Sz] = RCWA(layer,device,input_wave,wavelengthArray)
    arguments
        layer (:,1) struct
        device struct
        input_wave struct
        wavelengthArray (1,:) {mustBeNumeric,mustBeReal}
    end
    
c = onCleanup(@() progressBar());

if device.checkConverg
    harmArray = 1:2:device.Hmax;
else
    harmArray = device.Hmax;
end

issurf = cellfun(@(x) isa(x,"Surface"),{layer.input});
laynum = find(issurf,1);
if laynum
    k = [layer(issurf).input];
    k2 = [k.surfsize];
    k3 = [k.surfres];
    if ~all(k2 == k2(1))&&~all(k3 == k3(1))
        error("Mismatching layer dimensions")
    end
    
    if device.useSurfaceSize
        device.size_x       = layer(laynum).input.surfsize;
        device.size_y       = layer(laynum).input.surfsize;
        device.res_x        = layer(laynum).input.surfres;
        device.res_y        = layer(laynum).input.surfres;   
    end
elseif device.useSurfaceSize
    warning("No rough layers added, using manually entered dimensions")
end

lenHarmArray = length(harmArray);
lenWavelengthArray = length(wavelengthArray);

progressTick = progressBar(lenWavelengthArray*lenHarmArray);


for iHarm = 1:lenHarmArray
    
    device.P            = harmArray(iHarm);
    device.Q            = harmArray(iHarm);
    device.num_H        = device.P * device.Q;
    
    device              = truncation(device);
    input_wave          = get_sinc(input_wave,device);
    layer               = build_layerstack(layer,device);
    
    Sz = zeros(device.num_H,~device.calcAllRough*numel(layer)+device.calcAllRough*numel([layer.L])+2,lenWavelengthArray);

    parfor iWavelength = 1:lenWavelengthArray      
        
        Sz(:,:,iWavelength) = RCWA_transmittance(layer,device,input_wave,iWavelength);

        progressTick();
    end

R(:,iHarm) = squeeze(sum(Sz(:,end-1,:),1));

end

if device.checkConverg
    figure
    plot(wavelengthArray,R)
    %xticklabels(string(wavelengthArray))
    xlabel("Wavelength (nm)")
    ylabel("Reflectance (a.u.)")
    legend(string(harmArray.^2)+' Harmonics (excl. trunc)','location','eastoutside')
end

end


function progressOut = progressBar(varargin)
persistent iters wb tocArray;
if nargin==1
    i=varargin{1};

    D = parallel.pool.DataQueue;
    afterEach(D, @updateWaitbar);
    progressOut = @progressTick;
elseif nargin==0
    close(wb);
    delete(wb);
    toc
    clear wb iters tocArray
end

    function updateWaitbar(~)
        if isempty(wb)
            wb = waitbar(0,'Initializing...');
            iters = 1;
            tic;
            tocArray = 0;
        else
            tocArray(end+1) = toc;

            iterRemaining = i-iters;

            if iters > 10
                t = tocArray(end-10:end);
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

    function progressTick()
        send(D, []);
    end
end