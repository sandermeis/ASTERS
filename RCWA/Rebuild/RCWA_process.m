function RCWA_process(folderName)

p = load("results/"+folderName+"/param.mat","param");
param=p.param;

%algaasthickness = 1:1:10;
oxidethickness  = 180:40:220;
%[X,Y] = ndgrid(algaasthickness,oxidethickness);
%fulljsc = zeros(size(oxidethickness));

for i = 1:numel(param)
A = load("results/"+folderName+"/sim"+string(i)+".mat");
%fulljsc(i) = A.fom(1);
% SS = A.Sz;
% SS2 = squeeze(sum(SS,1));
% SS3(i,:) = SS2(4,:);
%plot(param(1).wavelengthArray,SS2(4,:),'LineWidth',2)
%hold on
%if numel(param)<11
    %warning("Too many simulations, not displaying")
    layer = fill_layer(param(i), "results/" + folderName);
    %displayDiscretized(layer(1).geometry.eps_struc, 2)
    try
    RCWA_plot(param(i), A.Sz, layer, i)
    catch
    end
%end
end
% plot(5:2:15,mean(abs(diff(SS3,1)),2)-0.02,'LineWidth',2)
% hold on
% bb=load('bb.mat');
% plot(9:2:19,bb.bb,'LineWidth',2)
% xlabel("Harmonics (#)","FontSize", 16, "FontWeight", 'bold')
% ylabel("Mean reflectivity difference (a.u.)","FontSize", 16, "FontWeight", 'bold')
% title("Convergence of two wedges of different sizes","FontSize", 16, "FontWeight", 'bold')
% legend(["2500 nm wide","7500 nm wide"])
%legend(["3","5","7","9","11","13","15"])

% figure('Color','w');
% plot(oxidethickness,fulljsc)
% 
% title("Jsc per simulation", "FontSize", 18, "FontWeight", 'bold')
% %xlabel("Algaas", "FontSize", 16, "FontWeight", 'bold')
% xlabel("Oxide", "FontSize", 16, "FontWeight", 'bold')
% ylabel("Jsc (mA/cm^2)", "FontSize", 16, "FontWeight", 'bold')


% h = RCWA_plot(fill_layer(param(3)),param(3),param(3).wavelengthArray, Sz)
%plot([param.oxidethickness],fulljsc)
%plot([param.res],fulljsc)
% 
% p = struct('p1', {param.p1}, 'p2', {param.p2},'p3', {param.p3},'p4', {param.p4});
% n = 176;
% %[layer] = fill_layer(p, n, [param.lay])
% %layer(8).input.plot
% %%
% % filter for param p1
% % average rest
% plist = [p.p3];
% [C,IA,IC] = unique(plist);
% avg = zeros(1,numel(C));
% for i=1:numel(C)
%     f = find(IC,i);
%     avg(i) = mean(fulljsc(f));
% end
% figure
% plot(C,avg)
%%
end



function disp_permittivity()
        % this can happen after everything is imported
        if options.dispPermeabilityFig
            t = tiledlayout(1,2);
            title(t,'Material properties of layer stack')
            nexttile
            hold on
            for i=1:length(ip)
                plot(wavelengthArray,ip(:,1),"LineWidth",2)
            end
            xlabel('Wavelength (nm)')
            ylabel('Refractive index')
            xlim([lab1,lab2])

            nexttile
            hold on
            for i=1:length(eps_lab)
                plot(wavelengthArray,ip(:,2),"LineWidth",2)
            end
            xlabel('Wavelength (nm)')
            ylabel('Extinction coefficient')
            xlim([lab1,lab2])
            legend(string({layer.material}),'location','eastoutside');
        end
end



function plotLayer(layer,i)
    %plot after processing is done
    if layer(i).input~=0
        figure % maybe something with reverse
        if layer(i).reverse
            mesh(sum(layer(i).geometry.eps_struc .* reshape(layer(i).L,1,1,[]),3))
            set(gca, 'zdir', 'reverse')
            zt = get(gca, 'ZTick');
            set(gca, 'ZTickLabel', fliplr(zt))
        else
            mesh(sum(layer(i).geometry.eps_struc .* reshape(layer(i).L,1,1,[]),3))
        end
        title("Layer "+i+": "+layer(i).material)
        zlabel("-Z")
    else
        warning('Uniform layer selected, cannot display this structure')
    end
end