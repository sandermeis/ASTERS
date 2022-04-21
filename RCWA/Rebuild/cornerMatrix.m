function C = cornerMatrix(A, dimvec, isParallel)
arguments
    A {isnumeric}
    dimvec {iscell} = {}
    isParallel {iscell} = {}
end
%%
% % select more dimensions to vary together
%dimvec = {};
%dimvec{1}=[2,3];
%diagdimvec{1}=[1,2,3];
% % %dimvec{2}=[2,3];
%dimvec{3}=[1,3];
%A = ones(3,3,4);

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
        if ~anymissing(dimvec{i})&&~isempty(dimvec{i})
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

