function folderName = RCWA(param, options)
arguments
    param
    options
end

dt = datestr(datetime,'dd_mm_yy_HH_MM_SS');
if ~isempty(options.simulationName)
    folderName = options.simulationName + "_" + dt;
else
    folderName = "sim_" + dt;
end

c = onCleanup(@() progressBar());


numRuns = numel(param);

fig = uifigure;
numSimWarning = uiconfirm(fig, sprintf("About to do %d simulations", numRuns), "title warning");
numCoreWarning = '';
if options.parallel
    numCores = feature('numCores');
    if numRuns <= numCores
        title = sprintf("%d simulations <= %d cores",numRuns,numCores);
        msg = "Non parallel might be slightly faster";
        numCoreWarning = uiconfirm(fig, msg, title, ...
            'Options', {'Parallel anyway', 'Non parallel', 'Cancel'}, ...
            'DefaultOption', 2, 'CancelOption', 3);
        if numCoreWarning == "Non parallel"
            options.parallel = 0;
        end
    end
end
if numSimWarning == "OK" && ~(numCoreWarning == "Cancel")
    close(fig)

    offlinePathName = "results/";
    mkdir(offlinePathName, folderName)
    offlinePathName = offlinePathName + folderName;
    save(offlinePathName + "/param.mat", "param")
    copyfile("input/input.txt", offlinePathName)
    copyfile("input/layers.xlsx", offlinePathName)
            % following line needs update, for loop over all
            % surfacefilenames
    copyfile("input/"+param(1).surfaceFile+".m", offlinePathName)

        if options.desktop
            onlinePathName = '/run/user/1000/gvfs/smb-share:server=amsbackup-srv.science.ru.nl,share=amsbackup/Students/Sander/results/';
        else
            onlinePathName = '/Volumes/amsbackup/Students/Sander/results/';
        end

    if options.onlinesave
        try
            mkdir(onlinePathName, folderName)
            onlinePathName = onlinePathName + folderName;
            copyfile("input/input.txt", offlinePathName)
            copyfile("input/layers.xlsx", offlinePathName)
            copyfile("input/"+param(1).surfaceFile+".m", offlinePathName)
            % following line needs update, for loop over all
            % surfacefilenames
            copyfile(offlinePathName + "/param.mat", onlinePathName)
        catch
            warning("Unable to create online directory")
            options.onlinesave = false;
        end
    end

    % REDO THIS, SKIPPING FOR NOW
    %     % check if surfaces same dimensions
    %     issurf = cellfun(@(x) isa(x,"Surface"),{layer.input});
    %     laynum = find(issurf,1);
    %     if laynum
    %         k = [layer(issurf).input];
    %         k2 = [k.surfsize];
    %         k3 = [k.surfres];
    %         if ~all(k2 == k2(1))&&~all(k3 == k3(1))
    %             error("Mismatching layer dimensions")
    %         end
    %
    %         if param.useSurfaceSize
    %             param.size_x       = layer(laynum).input.surfsize;
    %             param.size_y       = layer(laynum).input.surfsize;
    %             param.res_x        = layer(laynum).input.surfres;
    %             param.res_y        = layer(laynum).input.surfres;
    %         end
    %     elseif param.useSurfaceSize
    %         warning("No rough layers added, using manually entered dimensions")
    %     end

    [~] = progressBar(false, numel([param.wavelengthArray]));

    if options.parallel
        progressTick = progressBar(options.parallel);
        parfor n = 1:numRuns
            layer = fill_layer(param(n));
            [Sz, fields, Kx, Ky, Kz] = RCWA_transmittance(layer, param(n), progressTick);
            fom = Jsc(squeeze(sum(Sz,1)),param(n).wavelengthArray);

            fileName = "results/" + folderName + "/sim" + n + ".mat";

            if options.onlinesave
                parsave(fileName, onlinePathName, Sz, fields, Kx, Ky, Kz, fom, n)
            else
                parsave(fileName, [], Sz, fields, Kx, Ky, Kz, fom, n)
            end
        end
    else
        progressTick = progressBar(options.parallel);
        for n = 1:numRuns
            layer = fill_layer(param(n));
            [Sz, fields, Kx, Ky, Kz] = RCWA_transmittance(layer, param(n), progressTick);
            fom = Jsc(squeeze(sum(Sz, 1)), param(n).wavelengthArray);

            fileName = "results/" + folderName + "/sim" + n + ".mat";

            if options.onlinesave
                parsave(fileName, onlinePathName, Sz, fields, Kx, Ky, Kz, fom, n)
            else
                parsave(fileName, [], Sz, fields, Kx, Ky, Kz, fom, n)
            end
        end
    end
else
    close(fig)
end

end


function progressOut = progressBar(varargin)
persistent iters wb tocArray maxIter;
if nargin==2
    isPar = varargin{1};
    maxIter = varargin{2};

    wb = waitbar(0,'Working on first iteration...');
    iters = 1;
    tic;
    tocArray = 0;

    if isPar
        D = parallel.pool.DataQueue;
        afterEach(D, @updateWaitbar);
        progressOut = @progressTick;
    else
        progressOut = @updateWaitbar;
    end
elseif nargin==1
    isPar=varargin{1};
    if isPar
        D = parallel.pool.DataQueue;
        afterEach(D, @updateWaitbar);
        progressOut = @progressTick;
    else
        progressOut = @updateWaitbar;
    end
elseif nargin==0
    close(wb);
    delete(wb);
    toc
    clear wb iters tocArray
end

    function updateWaitbar(~)
        tocArray(end+1) = toc;

        iterRemaining = maxIter-iters;

        %             if iters > 10
        %                 t = tocArray(end-10:end);
        %             else
        t = tocArray;
        %             end

        timeLeft = string(seconds(iterRemaining*mean(diff(t))),'hh:mm:ss');

        waitBarDuration = iters/maxIter;
        waitBarString = {'Iteration ' + string(iters) + '/' + string(maxIter),...
            'Estimated time remaining: ' + timeLeft};

        waitbar(waitBarDuration,wb,waitBarString)

        iters = iters + 1;
    end

    function progressTick()
        send(D, []);
    end
end

function parsave(fileName, onlinePathName, Sz, fields, Kx, Ky, Kz, fom, n)
% savefile = varargin{1}; % first input argument
% for i=2:nargin
%     savevar.(inputname(i)) = varargin{i}; % other input arguments
% end
warning("Saving to file")
save(fileName, 'Sz', 'fields', 'Kx', 'Ky', 'Kz', 'fom', 'n')
if ~isempty(onlinePathName)
    try
        copyfile(onlinePathName, fileName)
    catch
        warning("Unable to save file online")
    end
end
end