function RCWA_process(folderName, whichsims)
arguments
folderName (1,:) {isstring}
whichsims = []
end

for j = 1:numel(folderName)

    % Loads parameter set
    p = load("results/" + folderName(j) + "/param.mat", "param");
    param = p.param;

    % Loop through parameter set
    for i = 1:numel(param)
        if isempty(whichsims) || ismember(i, whichsims)
            A = load("results/" + folderName(j) + "/sim" + string(i) + ".mat");
            layer = fill_layer(param(i), "results/" + folderName(j));
    
    
            % plot graphs of jsc per simulation, auto detect variables,
            % harmonics = convergence plot

            % figure('Color','w');
            % 
            % plot(szz,SS3)
            % title("Jsc per simulation", "FontSize", 18, "FontWeight", 'bold')
            % xlabel("size", "FontSize", 16, "FontWeight", 'bold')
            % ylabel("Jsc (mA/cm^2)", "FontSize", 16, "FontWeight", 'bold')
            % legend
    
            % displayDiscretized(layer(1).geometry.eps_struc, 2)

            RCWA_plot(param(i), A.Sz, layer, i, "Flat sim " + string(i), 4)
        end
    end

    %show_fields

    %plotlayer
    %disp permittivity

end
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