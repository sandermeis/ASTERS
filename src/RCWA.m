function folderName = RCWA(param, sim_options)
arguments
    param
    sim_options.simulationName = 'Simulation';
    sim_options.parallel = false;
    sim_options.onlinesave = false;
    sim_options.onlinePathName = '/run/user/1000/gvfs/smb-share:server=amsbackup-srv.science.ru.nl,share=amsbackup/Students/Sander/results/';
    sim_options.mode {mustBeMember(sim_options.mode,["ET", "S"])} = "ET"
end

% Set simulation name based on date and time
dt = datestr(datetime,'dd_mm_yy_HH_MM_SS');
if ~isempty(sim_options.simulationName)
    folderName = sim_options.simulationName + "_" + dt;
else
    folderName = "sim_" + dt;
end

c = onCleanup(@() progressBar());

mode = sim_options.mode;
numRuns = numel(param);

% Confirmation box
fig = uifigure;
numSimWarning = uiconfirm(fig, sprintf("About to do %d simulations", numRuns), "Start simulations");
numCoreWarning = '';
if sim_options.parallel
    numCores = feature('numCores');
    if numRuns <= numCores
        title = sprintf("%d simulations <= %d cores",numRuns,numCores);
        msg = "Non parallel might be slightly faster";
        numCoreWarning = uiconfirm(fig, msg, title, ...
            'Options', {'Parallel anyway', 'Non parallel', 'Cancel'}, ...
            'DefaultOption', 2, 'CancelOption', 3);
        if numCoreWarning == "Non parallel"
            sim_options.parallel = 0;
        end
    end
end

% Proceed
if numSimWarning == "OK" && ~(numCoreWarning == "Cancel")
    close(fig)
    
    % Create path in results folder
    offlinePathName = "results/";
    mkdir(offlinePathName, folderName)

    offlinePathName = offlinePathName + folderName;

    % Save parameter file, input, layers and surface file(s)
    save(offlinePathName + "/param.mat", "param")
    copyfile("input/input.txt", offlinePathName)
    copyfile("input/layers.xlsx", offlinePathName)
    surfaceFiles = [param.surfaceFile];
    for sf = 1:numel(surfaceFiles)
        copyfile("input/" + surfaceFiles(sf) + ".m", offlinePathName)
    end

    % If enabled, copy these files to online location
    if sim_options.onlinesave
        try
            mkdir(sim_options.onlinePathName, folderName)
            onlinePathName = sim_options.onlinePathName + folderName;

            % Save parameter files in online folder
            copyfile(offlinePathName + "/param.mat", onlinePathName)
            copyfile("input/input.txt", onlinePathName)
            copyfile("input/layers.xlsx", onlinePathName)

            surfaceFiles = [param.surfaceFile];
            for sf = 1:numel(surfaceFiles)
                copyfile("input/" + surfaceFiles(sf) + ".m", onlinePathName)
            end

        catch
            warning("Unable to create online directory")
            sim_options.onlinesave = false;
        end
    end

    % Initialize progress bar
    [~] = progressBar(false, numel([param.wavelengthArray]));
    
    % Simulation in parallel
    if sim_options.parallel
        progressTick = progressBar(sim_options.parallel);
        parfor n = 1:numRuns
            layer = fill_layer(param(n));
            [Sz, fields, Kx, Ky, Kz] = run_RCWA(layer, param(n), mode, progressTick);

            fileName = "results/" + folderName + "/sim" + n + ".mat";

            if sim_options.onlinesave
                parsave(fileName, onlinePathName, Sz, fields, Kx, Ky, Kz, n)
            else
                parsave(fileName, [], Sz, fields, Kx, Ky, Kz, n)
            end
        end
    % Simulation in series
    else
        progressTick = progressBar(sim_options.parallel);
        % Loop over parameter set
        for n = 1:numRuns

            % Create layer struct using parameters
            layer = fill_layer(param(n));

            % Run simulation
            [Sz, fields, Kx, Ky, Kz] = run_RCWA(layer, param(n), mode, progressTick);

            fileName = "results/" + folderName + "/sim" + n + ".mat";

            % Save these parameters as .mat in previously created folder
            if sim_options.onlinesave
                parsave(fileName, onlinePathName, Sz, fields, Kx, Ky, Kz, n)
            else
                parsave(fileName, [], Sz, fields, Kx, Ky, Kz, n)
            end
        end
    end
else
    close(fig)
end

end


% Progress bar which also works in parallel simulations
function progressOut = progressBar(varargin)
persistent iters wb tocArray maxIter;
if nargin == 2
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
elseif nargin == 1
    isPar=varargin{1};
    if isPar
        D = parallel.pool.DataQueue;
        afterEach(D, @updateWaitbar);
        progressOut = @progressTick;
    else
        progressOut = @updateWaitbar;
    end
elseif nargin == 0
    close(wb);
    delete(wb);
    toc
    clear wb iters tocArray
end

    function updateWaitbar(~)
        tocArray(end + 1) = toc;

        iterRemaining = maxIter - iters;

        t = tocArray;

        timeLeft = string(seconds(iterRemaining * mean(diff(t))),'hh:mm:ss');

        waitBarDuration = iters / maxIter;
        waitBarString = {'Iteration ' + string(iters) + '/' + string(maxIter),...
            'Estimated time remaining: ' + timeLeft};

        waitbar(waitBarDuration,wb,waitBarString)

        iters = iters + 1;
    end

    function progressTick()
        send(D, []);
    end
end


% Save results during execution
function parsave(fileName, onlinePathName, Sz, fields, Kx, Ky, Kz, n)
warning("Saving to file")
save(fileName, 'Sz', 'fields', 'Kx', 'Ky', 'Kz', 'n')
if ~isempty(onlinePathName)
    try
        copyfile(onlinePathName, fileName)
    catch
        warning("Unable to save file online")
    end
end
end