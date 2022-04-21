function [param] = load_parameters()

cellpar = readcell("input/input.txt", 'CommentStyle', '%','Delimiter',{'=','+'},'LineEnding',{'\n',';'},'FileType','text','TextType','string','DurationType','text');
numParams = size(cellpar,1);

% Parse input file
isListCells = false(numParams,1);
isScalarListCell = false(numParams,1);

for i = 1:numParams
    if isstring(cellpar{i,2})

        % Remove "&" and count
        if contains(cellpar{i,2},"&")
            cellpar{i,2} = erase(cellpar{i,2},"&");
            isScalarListCell(i) = 1;
        end

        cellpar{i,2} = parseText(cellpar{i,2});

        if size(cellpar,2)>2 && ~ismissing(cellpar{i,3})
            cellpar{i,3} = string(strsplit(cellpar{i,3},{' ',','}));
        end

        % Cells with array
        isListCells(i) = length(cellpar{i,2})>1;
    end
end
%%



% Remove "\" and count
% if contains(cellpar{i,2},"\")
%     cellpar{i,2} = erase(cellpar{i,2},"\");
%     cellSteady(i) = 1;
% end

% Exclude lists that should be considered as a scalar parameter
isListCells = isListCells & ~isScalarListCell;

% Indices of the cells with array
k = find(isListCells);
% Cell with just the arrays
ListCells = cellpar(isListCells, :);

dim = numel(ListCells(:,2));

% Process extra parameters
if size(cellpar, 2)>2
    for j = 1:numel(ListCells(:, 3))
        % For non missing arrays
        if all(~ismissing(ListCells{j, 3}))
            d1 = [];
            isParallel = [];
            isRange = {};
            for i = 1:numel(ListCells{j,3}) % loop through combinations
                if contains(char(ListCells{j,3}(i)),"(")
                    ss = strsplit(ListCells{j,3}(i),"(");
                    ss(2) = "("+ss(2);
                    isRange{1} = {1};
                    isRange{end+1} = {parseText(ss(2))};
                    ListCells{j,3}(i) = ss(1);

                    if contains(char(ListCells{j,3}(i)),"\")
                        ListCells{j,3}(i) = erase(ListCells{j,3}(i),"\");
                        isParallel(1) = true;
                        isParallel(end+1) = true;
                    else
                        isParallel(1) = false;
                        isParallel(end+1) = false;
                    end
                else
                    if contains(char(ListCells{j,3}(i)),"\")
                        ListCells{j,3}(i) = erase(ListCells{j,3}(i),"\");
                        isParallel(1) = true;
                        isParallel(end+1) = true;
                    else
                        isParallel(1) = false;
                        isParallel(end+1) = false;
                    end
                end
                % Replace with vector
                d1(1) = j;
                paramNames = ListCells(:,1);
                paramLoc = find([paramNames{:}]==ListCells{j,3}(i));
                if ~isempty(paramLoc)
                    d1(end+1) = paramLoc;
                else
                    pname = string(ListCells{j,3}(i));
                    error("Unable to find parameter " + pname)
                end
            end
            ListCells{j,3} = d1;
            ListCells{j,4} = isParallel;
            ListCells{j,5} = isRange;
        end
    end
end

b = cell(1,dim);

if dim>0
    % Make grid with all lists
    [b{:}] = ndgrid(ListCells{:,2});
    
    if size(cellpar,2)>2
        cm = cornerMatrix(b{1},ListCells(cellfun(@(x) isnumeric(x), ListCells(:,3)),3), ListCells(cellfun(@(x) isnumeric(x), ListCells(:,4)),4));
    else % Only vertices
        cm = cornerMatrix(b{1},{});
    end

    for i=1:dim
        c{i} = b{i}(cm);
    end
    
    for n = 1:length(c{1,1})
        cp2 = cellpar;
        for j = 1:dim
            %if numeric
            if iscell(c{j}(n))
                cp2(k(j),2)=c{j}(n);
            else
                cp2(k(j),2)={c{j}(n)};
            end
        end

        param(n) = cell2struct(cp2(:,2), string(cp2(:,1)), 1);
    end
else
    param(1) = cell2struct(cellpar(:,2), string(cellpar(:,1)), 1);
end
% Do always

%%

[param.size_x]       = deal(param.size);
[param.size_y]       = deal(param.size);
[param.res_x]        = deal(param.res);
[param.res_y]        = deal(param.res);

test                    = cellfun(@(x) 2*pi./x,{param.wavelengthArray},'UniformOutput',false);
[param.k_0]          = deal(test{:});

[param.P]            = deal(param.Hmax);
[param.Q]            = deal(param.Hmax);
[param.num_H]        = deal([param.P] .* [param.Q]);
%%

param     = truncation(param);

% to not repeat look ups:
% convert list to string, place next to eachother, and find unique by row
[C_ref, ia_ref, ic_ref] = unique([[param.ref_medium]',string(cellfun(@(x) num2str(x),{param.wavelengthArray},'UniformOutput',false))'],'rows');

eps_lab_ref   = import_permittivities({C_ref(:,1)});%, {param(ia_ref).wavelengthArray});
wl_ref = {param(ia_ref).wavelengthArray};

% {1} for interoperability with fill_layer
for i=1:numel(eps_lab_ref{1})
    eps_lab_ref2{i}   = eps_lab_ref{1}{i}(wl_ref{i});
end

[param.eps_ref]     = deal(eps_lab_ref2{ic_ref});
throwaway           = cellfun(@(x) ones(size(x)),{param(:).wavelengthArray},'UniformOutput',false);
[param.mu_ref]      = deal(throwaway{:});

[C_trn, ia_trn, ic_trn] = unique([[param.trn_medium]',string(cellfun(@(x) num2str(x),{param.wavelengthArray},'UniformOutput',false))'],'rows');
eps_lab_trn   = import_permittivities({C_trn(:,1)});%, {param(ia_trn).wavelengthArray});
wl_trn = {param(ia_trn).wavelengthArray};

% {1} for interoperability with fill_layer
for i=1:numel(eps_lab_trn{1})
    eps_lab_trn2{i}   = eps_lab_trn{1}{i}(wl_trn{i});
end

[param.eps_trn]     = deal(eps_lab_trn2{ic_trn});
[param.mu_trn]      = deal(throwaway{:});


end




% a = Surface(512,10000);
% f_titan = Feature(round(512/2),5000,500,"Cone");
% f_giant = Feature(round(512/2.5),5000,450,"Cone");
% f_large = Feature(round(512/3.5),2500,400,"Cone");
% f_medium = Feature(round(512/5),1250,300,"Cone");
% f_small = Feature(round(512/10),625,100,"Cone");
%
% a.addRandomFeatures(f_titan,2,"PBC",true)
% a.addRandomFeatures(f_giant,5,"PBC",true)
% a.addRandomFeatures(f_large,15,"PBC",true)
% a.addRandomFeatures(f_medium,80,"PBC",true)
% a.addRandomFeatures(f_small,100,"PBC",true)
%
% a.placeFeatures("PBC",true, "mode", "merge");
% a.listFeatures()
% a.report()
% a.plot()
% a.resize(128,5000)
% a.plot()
% b4 = Feature("Feature1_2500nm.csv",2500);
%
% b5 = Feature("Feature1_2500nm.csv",2500);
% b5 = rotate(b5);
%a.addFeature(b,1,1);
% a.addFeature(b5,1,1)
% a.addFeature(b2,1,1)
% a.addFeature(b3,1,1)
%

% a.placeRandomFeatures(1,1000,"PBC",true, "mode", "merge")
%a.hscale(200)
%a.addRandomFeatures(b4,3,"PBC",true)
%a.addRandomFeatures(b5,2,"PBC",true)
% a.placeFeatures("PBC",true, "mode", "merge");
%a.placeRandomFeatures(2,3,"PBC",true, "mode", "merge")
%a.plot();

% % Constants
% param.version      = "3.0.3";
%
% % Parameters
% wavelengthArray     = 350:5:900;
% param.ref_medium	= "Air";
% param.trn_medium	= "Air";
% param.theta         = 0;
% param.phi           = 0;
% param.pTE           = 0.5;
% param.pTM           = 0.5;
% param.size_x       = 10000;
% param.size_y       = 10000;
% param.res_x        = 512;
% param.res_y        = 512;
% param.truncFactor  = 1;
% param.Hmax         = 9;
% param.tolerance    = 5e-3;
%
% % Booleans
% param.plot_permfig     = false;
% param.useSurfaceSize   = false;
% param.truncfig         = false;
% param.plotSurf         = false;


% param.truncateHarm     = false;

% param.calcAllRough     = false;

% param.reverse          = false;
% param.recalcRoughL     = false;
% param.optimRough       = false;

%options.dispTruncFig
%options.dispPermeabilityFig
%options.dispSurface
%options.checkConvergence + which dimension
%%

function output = parseText(input)
% If it contains ',', parse text as list, else as scalar
output = parseEval(input);

if contains(input, ',')&&isequal(output,input)
    b = strsplit(input,{' ',','});
    output = cell(size(b));
    % need to add safety check if numeric
    for j = 1:numel(b)
        output{j} = parseEval(b(j));
    end
else
    output = parseEval(input);
end


    function out = parseEval(input)
        if all(ismember(char(input), '0123456789+-.:*/pi()[],')) ...
                || input == "false" || input == "true"
            try
                out = eval(input);
                if isempty(out)
                    error("Unknown numeric input: " + string(input))
                end
            catch
                error("Unknown numeric input: " + string(input))
                out = input;
            end
        else
            out = input;
        end
    end

end


function param = truncation(param)
%%
for i=1:numel(param)
    M = -(param(i).P-1)/2:(param(i).P-1)/2;
    N = -(param(i).Q-1)/2:(param(i).Q-1)/2;
    [m,n] = meshgrid(M,N);

    if param(i).truncateHarm
        TMAP = abs(m/((param(i).P-1)/2)).^(2*param(i).truncFactor) + abs(n/((param(i).Q-1)/2)).^(2*param(i).truncFactor);
        TMAP(isnan(TMAP))=1;
        TMAP = (TMAP <= 1);
        if param(i).truncfig
            figure
            imagesc(TMAP);
        end
    else
        TMAP = ones(size(m));
    end

    % Extract Array Indices
    param(i).tr_ind = find(TMAP(:));
    param(i).num_H = length(param(i).tr_ind);
end


end