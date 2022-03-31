function [param,fulljsc] = RCWA_process(folderName)

p = load("results/"+folderName+"/param.mat","param");
param=p.param;

for i = 1:numel(param)
jsc = load("results/"+folderName+"/sim"+string(i)+".mat","fom");
fulljsc(i) = jsc.fom(4);

end

plot(1:numel(param),fulljsc)

p = struct('p1', {param.p1}, 'p2', {param.p2},'p3', {param.p3},'p4', {param.p4});
n = 176;
%[layer] = fill_layer(p, n, [param.lay])
%layer(8).input.plot
%%
% filter for param p1
% average rest
plist = [p.p3];
[C,IA,IC] = unique(plist);
avg = zeros(1,numel(C));
for i=1:numel(C)
    f = find(IC,i);
    avg(i) = mean(fulljsc(f));
end
figure
plot(C,avg)
%%
end

