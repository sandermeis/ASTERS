function [param] = load_parameters()

cellpar = readcell("./input/input.txt", 'CommentStyle', '%','Delimiter',{'=','+'},'LineEnding',{'\n',';'},'FileType','text','TextType','string','DurationType','text');
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
            isParallel = [false];
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
                        isParallel(end+1) = false;
                    end
                else
                    if contains(char(ListCells{j,3}(i)),"\")
                        ListCells{j,3}(i) = erase(ListCells{j,3}(i),"\");
                        isParallel(1) = true;
                        isParallel(end+1) = true;
                    else
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

nd_array = cell(1,dim);

if dim>0
    % Produce an N-dimensional grid with all variable combinations
    [nd_array{:}] = ndgrid(ListCells{:, 2});
    
    % cornerMatrix produces a logical array in N-dimensional space to
    % determine the requested variable combinations.

    % If multiple variables are codependent (faces or volumes in ND space)
    if size(cellpar,2) > 2
        cm = cornerMatrix(nd_array{1}, ListCells(cellfun(@(x) isnumeric(x), ListCells(:, 3)), 3), ListCells(cellfun(@(x) islogical(x), ListCells(:, 4)), 4));
    % If variables are independent (only vertices in ND space)
    else
        cm = cornerMatrix(nd_array{1}, {});
    end

    for i = 1:dim
        c{i} = nd_array{i}(cm);
    end
    
    for n = 1:length(c{1, 1})
        cp2 = cellpar;
        for j = 1:dim
            %if numeric
            if iscell(c{j}(n))
                cp2(k(j), 2) = c{j}(n);
            else
                cp2(k(j), 2) = {c{j}(n)};
            end
        end

        param(n) = cell2struct(cp2(:, 2), string(cp2(:, 1)), 1);
    end
else
    param(1) = cell2struct(cellpar(:, 2), string(cellpar(:, 1)), 1);
end

%%
for i = 1:dim
    for f = 1:numel(param)
        param(f).bk(i).name = ListCells{i, 1};
        aaa = [param.(param(f).bk(i).name)];
        for j = 1:numel(ListCells{i, 2})
            temp_list(j) = find(aaa==ListCells{i, 2}(j), 1);
        end
        param(f).bk(i).val = temp_list;
        param(f).bk(i).list = temp_list;
    end
end

%% Fill

[param.size_x]       = deal(param.size);
[param.size_y]       = deal(param.size);
[param.res_x]        = deal(param.res);
[param.res_y]        = deal(param.res);

test                 = cellfun(@(x) 2*pi./x,{param.wavelengthArray},'UniformOutput',false);
[param.k_0]          = deal(test{:});

[param.P]            = deal(param.Hmax);
[param.Q]            = deal(param.Hmax);
[param.num_H]        = deal([param.P] .* [param.Q]);

%% Truncate
param = truncation(param);

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
throwaway           = cellfun(@(x) ones(size(x)), {param(:).wavelengthArray}, 'UniformOutput', false);
[param.mu_ref]      = deal(throwaway{:});

[C_trn, ia_trn, ic_trn] = unique([[param.trn_medium]', string(cellfun(@(x) num2str(x), {param.wavelengthArray},'UniformOutput',false))'],'rows');
eps_lab_trn   = import_permittivities({C_trn(:,1)});%, {param(ia_trn).wavelengthArray});
wl_trn = {param(ia_trn).wavelengthArray};

% {1} for interoperability with fill_layer
for i=1:numel(eps_lab_trn{1})
    eps_lab_trn2{i}   = eps_lab_trn{1}{i}(wl_trn{i});
end

[param.eps_trn]     = deal(eps_lab_trn2{ic_trn});
[param.mu_trn]      = deal(throwaway{:});


end

function output = parseText(input)
% If it contains ',', parse text as list, else as scalar
output = parseEval(input);

if contains(input, ',') && isequal(output,input)
    b = strsplit(input, {' ',','});
    output = cell(size(b));
    % need to add safety check if numeric
    for j = 1:numel(b)
        output{j} = parseEval(b(j));
    end
else
    output = parseEval(input);
end


    function out = parseEval(input)
        if all(ismember(char(input), ' 0123456789+-.:*/pi()[],')) ...
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
for i = 1:numel(param)
    M = - (param(i).P - 1) / 2:(param(i).P - 1) / 2;
    N = - (param(i).Q - 1) / 2:(param(i).Q - 1) / 2;
    [m, n] = meshgrid(M, N);

    if param(i).truncateHarm
        TMAP = abs(m/((param(i).P-1)/2)).^(2*param(i).truncFactor) + abs(n/((param(i).Q-1)/2)).^(2*param(i).truncFactor);
        TMAP(isnan(TMAP))=1;
        TMAP = (TMAP <= 1);
    else
        TMAP = ones(size(m));
    end

    % Extract Array Indices
    param(i).tr_ind = find(TMAP(:));
    param(i).num_H = length(param(i).tr_ind);
end


end

function C = cornerMatrix(A, dimvec, isParallel)
arguments
    A {isnumeric}
    dimvec {iscell} = {}
    isParallel {iscell} = {}
end
%%
% % select more dimensions to vary together
% dimvec{1}=[1,2,3];
% A = ones(3,3,3);
% isParallel{1} = [1,1,0];

% Space dimensions
sz = size(A);
C = zeros(sz);

% Initialize cell
id = cell(1, numel(sz));

% Indices of each dimension 
[id{:}] = ind2sub(sz, find(ones(sz)));

% Matrix with dimension indices
cid = cell2mat(id);

% Vertices

% Dimension index = 1; so anywhere any dim is 1, all "Faces"
zid = (cid==1);

% if ~isempty(isParallel)
%     for i=1:numel(isParallel)
%         if ~anymissing(isParallel{i}) && ~isempty(isParallel{i})
%             temp = zid;
%             temp2 = cid;
%             temp(:, isParallel{i}) = 1;
%             % Definitely include input dimensions
%             a = temp2(:, isParallel{i});
%             a2 = all(diff(a,1,2)==0,2);
%             % All indices 1 except those in dimvec, so make those 1
%             bid = all([temp==1,a2],2);%all([any(zid(:, dimvec{i})==0, 2), all(temp, 2)], 2);
%             %B = all(cid,2);
%             C(bid) = 1;
%         end
%     end
% end

% Find only the "Faces" you want to vary together
if ~isempty(dimvec)
    for i=1:numel(dimvec)
        if ~any(ismissing(dimvec{i}))&&~isempty(dimvec{i})
            temp = zid;
            temp2 = cid;
            whichParallel{i} = dimvec{i}(logical(isParallel{i}));
            a = temp2(:, whichParallel{i});
            a2 = all(diff(a,1,2)==0,2);
            % Definitely include input dimensions
            temp(:, dimvec{i}) = 1;

            % All indices 1 except those in dimvec, so make those 1
            bid = all([temp==1,a2],2);%all([any(zid(:, dimvec{i})==0, 2), all(temp, 2)], 2);
            %B = all(cid,2);
            C(bid) = 1;
        end
    end
end

% Of this set, take just the coordinates where all are one, except one, the "Vertices"
rid = (sum(zid, 2) >= (numel(sz)-1));
% Exclude parallel dimensions, where any

if ~isempty(isParallel)
    upar = unique([whichParallel{:}]);
    corrected_rid = rid & all(zid(:,upar),2);
    % Write vertices to output
    C(corrected_rid) = 1;
else
    C(rid) = 1;
end

% Turn the array into logical indexing
C = logical(C);
%%
end
