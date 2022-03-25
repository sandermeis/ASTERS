function C = cornerMatrix(A,dimvec)
arguments
    A {isnumeric}
    dimvec {iscell} = {}
end
%%
% % select more dimensions to vary together
% dimvec{1}=[1,2,3,4];
% % %dimvec{2}=[2,3];
% % %dimvec{3}=[1,3];
% A=ones(3,3,3,2);

sz=size(A);
C=zeros(sz);
id=cell(1,numel(sz));
[id{:}] = ind2sub(sz,find(ones(sz)));
cid=cell2mat(id);

%corners
zid=(cid==1);
rid=(sum(zid,2)>=(numel(sz)-1));
C(rid)=1;

%faces
if ~isempty(dimvec)
    for i=1:numel(dimvec)
        if ~ismissing(dimvec{i})
    temp = zid;
    temp(:,dimvec{i})=1;
    bid=all([any(zid(:,dimvec{i})==0,2),all(temp,2)],2);
    %B = all(cid,2);
    C(bid)=1;
        end
    end
end
C = logical(C);
%%
end

