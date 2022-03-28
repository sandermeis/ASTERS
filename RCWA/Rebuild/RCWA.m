function Sz = RCWA(layer_structure, param, options)
arguments
    layer_structure (1,:) struct
    param struct
    options
end

%c = onCleanup(@() progressBar());

if options.save
    folderName = "sim_" + datestr(datetime,'dd_mm_yy_HH_MM_SS');
% onlinepath='schijf/sander/results';
    mkdir("results",folderName)
    save("results/"+folderName+"/param.mat","param","layer_structure")
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

numRuns = numel(param);

%progressTick = progressBar(numRuns);

for n = 1:numRuns

    Sz = RCWA_transmittance(layer_structure(param(n).lay).layer_set(param(n).pset).layer, param(n));
    fom = Jsc(squeeze(sum(Sz,1)),param(n).wavelengthArray);

    if options.save
        fileName = "results/"+folderName+"/sim"+n+".mat";
        parsave(fileName,Sz,fom,n)
    end

    %progressTick();
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

function parsave(varargin)
savefile = varargin{1}; % first input argument
for i=2:nargin
    savevar.(inputname(i)) = varargin{i}; % other input arguments
end
save(savefile,'-struct','savevar')
end