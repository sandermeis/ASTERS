function [param] = load_parameters()

cellpar = readcell("input.txt", 'CommentStyle', '%','Delimiter',{'=','+','-'},'LineEnding',{'\n',';'},'FileType','text','TextType','char','DurationType','text');
numParams = size(cellpar,1);

for i = 1:numParams
    if ~any(ismember(cellpar{i,2}, ','))
        if ~isnumeric(cellpar{i,2})
        if all(ismember(char(cellpar{i,2}), '0123456789+-.:*/pi()[]'))||string(cellpar{i,2})=="false"||string(cellpar{i,2})=="true"
            try
            cellpar{i,2} = eval(cellpar{i,2});
            catch
            end
        else
            cellpar{i,2} = string(cellpar{i,2});
        end
        end
    else
        b = string(strsplit(cellpar{i,2},{' ',','}));
        b2=cell(size(b));
        % need to add safety check if numeric
        for j = 1:numel(b)
            try
            if all(ismember(char(b(j)), '0123456789+-.:[] '))||string(b(j))=="false"||string(b(j))=="true"
                b2{j} = eval(b(j));
            else
                b2{j} = b(j);
            end
            catch
            end
        end
        cellpar{i,2} = b2;
    end
    try
    cellpar{i,3} = string(strsplit(cellpar{i,3},{' ',','}));
    catch
    end
end

%cells with array, except with single c as parameter
cellsNumAndList = false(numParams,1);
if size(cellpar,2)>2
    for i = 1:numParams
        cellsNumAndList(i) = length(cellpar{i,2})>1 && any(~strcmp(cellpar{i,3},'c'));%isnumeric(cellpar{i,2})&&
    end
else
    for i = 1:numParams
        cellsNumAndList(i) = length(cellpar{i,2})>1;
    end
end
%indices cells with array
k = find(cellsNumAndList);
%cell with just the arrays
cc = cellpar(cellsNumAndList,:);
if size(cellpar,2)>2
for j=1:numel(cc(:,3)) %number of arrays
    if all(~ismissing(cc{j,3}))&&all(~strcmp(cc{j,3},'c')) % for non missing arrays
        d=[];
        for i = 1:numel(cc{j,3}) % loop through combinations
            d(1) = j;
            d(1+i) = find(cellfun(@(x) ~isempty(x),strfind(cc(:,1),cc{j,3}(i))));
        end
        cc{j,3} = d;
    end
end
end
dim = numel(cc(:,2));
b = cell(1,dim);
[b{:}] = ndgrid(cc{:,2});

%if params dont vary at same time
if size(cellpar,2)>2
    cm = cornerMatrix(b{1},cc(cellfun(@(x) isnumeric(x),cc(:,3)),3));
else
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
        %if not numeric
        %cp2(k(j),2)={cellpar(c{j}(n))};
    end
    param(n) = cell2struct(cp2(:,2),cp2(:,1),1);
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

eps_lab_ref   = import_permittivities(cellstr(C_ref(:,1)));%, {param(ia_ref).wavelengthArray});
wl_ref = {param(ia_ref).wavelengthArray};
for i=1:numel(eps_lab_ref)
    eps_lab_ref2{i}   = eps_lab_ref{i}(wl_ref{i});
end

[param.eps_ref]     = deal(eps_lab_ref2{ic_ref});
throwaway           = cellfun(@(x) ones(size(x)),{param(:).wavelengthArray},'UniformOutput',false);
[param.mu_ref]      = deal(throwaway{:});

[C_trn, ia_trn, ic_trn] = unique([[param.trn_medium]',string(cellfun(@(x) num2str(x),{param.wavelengthArray},'UniformOutput',false))'],'rows');
eps_lab_trn   = import_permittivities(cellstr(C_trn(:,1)));%, {param(ia_trn).wavelengthArray});
wl_trn = {param(ia_trn).wavelengthArray};
for i=1:numel(eps_lab_trn)
    eps_lab_trn2{i}   = eps_lab_trn{i}(wl_trn{i});
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