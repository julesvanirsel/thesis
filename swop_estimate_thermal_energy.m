do_plot = true;
save_plot = true;
happy = true;

direc = fullfile('data', 'paper2');
save_direc = fullfile(direc, 'dasc_data', 'thermal_energies');
event_fn = fullfile(direc, 'event_data.txt');
lines = readlines(event_fn)';
Qfs = 0.9:-0.1:0.1;
redfs = 0.1:0.1:0.4;
thermal_energies_fitted = nan(length(Qfs), length(redfs), length(lines)-2);
thermal_energies_maxima = thermal_energies_fitted;
thermal_energies_gofs = thermal_energies_fitted;
thermal_energies_error = thermal_energies_fitted;

choice = [4, 3];
if happy
    iqs = choice(1);
    irs = choice(2);
else
    iqs = 1:length(Qfs);
    irs = 1:length(redfs);
end

%%
for iq = iqs
    disp(iq)
    for ir = irs
        Q_filter_quantile = Qfs(iq);
        red_filter_quantile = redfs(ir);
        for il = 1:length(lines)-2
            line = lines(il+1);
            data = strsplit(line);
            time = datetime(data(2), 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
            time.Format = 'uuuuMMdd';
            h5fn = sprintf('INV_PKR_%s_%i.h5', time, second(time, 'secondofday'));
            h5fn = fullfile(direc, 'dasc_data', h5fn);
            glat = h5read(h5fn, '/Coordinates/Latitude');
            glon = h5read(h5fn, '/Coordinates/Longitude');
            Q = h5read(h5fn, '/Derived/Energy/Flux');
            E0 = h5read(h5fn, '/Derived/Energy/Characteristic');
            red = h5read(h5fn, '/Optical/Red');
            glat(isnan(E0)) = nan;
            glon(isnan(E0)) = nan;
            
            % filter of large Q and dim red
            E0_filtered = E0;
            E0_filtered(Q > quantile(Q(:), Q_filter_quantile)) = nan;
            E0_filtered(red < quantile(red(:), red_filter_quantile)) = nan;
            
            % determine energy histograms
            bin_edges = linspace(10, 3e3, 301);
            bin_means = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
            hist1 = histcounts(E0, bin_edges);
            hist2 = histcounts(E0_filtered, bin_edges);

            % determine thermal energies by maximum
            [m1, id1] = max(hist1);
            [m2, id2] = max(hist2);
            temp_max_unfiltered = bin_means(id1);
            temp_max_filtered = bin_means(id2);
            
            % determine thermal energies through gaussian fits
            start1 = [exp(1)*m1 / temp_max_unfiltered^2, temp_max_unfiltered];
            start2 = [exp(1)*m2 / temp_max_filtered^2, temp_max_filtered];
            [f1, gof1] = fit(bin_means', hist1', 'a*x^2*exp(-(x/b)^2)', 'start', start1);
            [f2, gof2] = fit(bin_means', hist2', 'a*x^2*exp(-(x/b)^2)', 'start', start2);
            temp_fit_unfiltered = f1.b;
            temp_fit_filtered = f2.b;

            % record data
            thermal_energies_fitted(iq, ir, il) = temp_fit_filtered;
            thermal_energies_maxima(iq, ir, il) = temp_max_filtered;
            thermal_energies_gofs(iq, ir, il) = gof2.adjrsquare;
            f2conf95 = confint(f2, 0.95);
            f2db = f2conf95(2,2) - temp_fit_filtered;
            thermal_energies_error(iq, ir, il) = f2db;
            
            if do_plot
                close all
                figure('Position', [10, 60, 1000, 800], 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 7, 10])
                tlo = tiledlayout(2, 2, 'TileSpacing','tight', 'Padding','tight');
                title(tlo, sprintf('Q filtered out: top %0i%% quantile, Red filtered out: lower %0i%% quantile, Calculated source region thermal energy: %i±%i eV', ...
                    round(100*(1-Q_filter_quantile)), round(100*red_filter_quantile), round(temp_fit_filtered, -1), round(f2db, -1)))

                nexttile
                pcolor(glon, glat, Q)
                shading flat
                colorbar
                clim([0, 20])
                xlabel('Geodetic longitude (°)')
                ylabel('Geodetic latitude (°)')
                title('Q (mW/m^2)')

                nexttile
                pcolor(glon, glat, E0)
                shading flat
                colorbar
                clim([0, 3e3])
                xlabel('Geodetic longitude (°)')
                ylabel('Geodetic latitude (°)')
                title('Unaccelerated Maxwellian E0 (eV)')

                nexttile
                pcolor(glon, glat, E0_filtered)
                shading flat
                colorbar
                clim([0, 3e3])
                xlabel('Geodetic longitude (°)')
                ylabel('Geodetic latitude (°)')
                title('Filtered Unaccelerated Maxwellian E0 (eV)')

                nexttile
                hold on
                histogram(E0, bin_edges, 'FaceColor', [1/2, 1/2, 1], 'EdgeColor', 'none');
                histogram(E0_filtered, bin_edges, 'FaceColor', [1, 1/2, 1/2], 'EdgeColor', 'none');
                plot(bin_means, f1(bin_means), 'b', 'LineWidth', 2)
                plot(bin_means, f2(bin_means), 'r', 'LineWidth', 2)
                plot([temp_fit_unfiltered, temp_fit_unfiltered], [0, f1(f1.b)], 'Color', [0, 0, 1], 'LineWidth', 1)
                plot([temp_max_unfiltered, temp_max_unfiltered], [0, m1], 'Color', [0, 0, 1], 'LineWidth', 1, 'LineStyle', '--')
                plot([temp_fit_filtered, temp_fit_filtered], [0, f2(f2.b)], 'Color', [1, 0, 0], 'LineWidth', 1)
                plot([temp_max_filtered, temp_max_filtered], [0, m2], 'Color', [1, 0, 0], 'LineWidth', 1, 'LineStyle', '--')
                xlim([0, 3e3]); ylim([0, max([hist1, hist2])*1.05])
                xlabel('Thermal energy, E0 (keV)')
                legend('Unfiltered', 'Filtered', ...
                    sprintf('Unfiltered fit (%i eV peak)', round(f1.b, -1)), ...
                    sprintf('Filtered fit (%i eV peak)', round(f2.b, -1)))
                
                if save_plot
                    filename = sprintf('%i%02i%02i_%i.png', ...
                        year(time), month(time), day(time), second(time, 'secondofday'));
                    filename = fullfile(save_direc, filename);
                    fprintf('Saving %s\n', filename)
                    exportgraphics(gcf, filename, 'Resolution', 600)
                    close all
                else
                    input('Press any key to continue ... ')
                end
                fprintf('Thermal energy = %i eV for %s\n', round(temp_fit_filtered, -1), time)
            else
                time.Format = "default";
                fprintf('Thermal energy = %i eV for %s\n', round(temp_fit_filtered, -1), time)
            end
        end
    end
end
%%

save('data\paper2\dasc_data\thermal_energies\filter_dependancies.mat', ...
    'thermal_energies_fitted', 'thermal_energies_maxima', ...
    'thermal_energies_error', 'thermal_energies_gofs')

%%
load('data\paper2\dasc_data\thermal_energies\filter_dependancies.mat')

close all
time.Format = 'uuuu-MM-dd';
qlim = 8;

figure('Position', [10, 60, 1500, 800])
tiledlayout(3,4)
for ie = 1:11
    nexttile
    hold on
    tmp = squeeze(thermal_energies_fitted(1:qlim,:,ie));
    imagesc(tmp)
    scatter(choice(2), choice(1), 10, 'r', 'filled')
    colorbar
    clim(thermal_energies_fitted(choice(1),choice(2),ie) + [-1,1]*200)
    xlim([0.5,0.5+size(tmp,2)]); ylim([0.5,0.5+size(tmp,1)])
    xlabel('Bottom red quantile removed (%) / 10')
    ylabel('Top Q quantile removed (%) / 10')
    title(sprintf('%s, %i eV spread', time, round(range(tmp(:)))))
end
exportgraphics(gcf, fullfile(save_direc, 'thermal_energies_fitted.png'), 'Resolution', 600)

figure('Position', [10, 60, 1500, 800])
tiledlayout(3,4)
for ie = 1:11
    nexttile
    hold on
    tmp = squeeze(thermal_energies_gofs(1:qlim,:,ie));
    imagesc(tmp)
    scatter(choice(2), choice(1), 10, 'r', 'filled')
    colorbar
    clim([0,1])
    xlim([0.5,0.5+size(tmp,2)]); ylim([0.5,0.5+size(tmp,1)])
    xlabel('Bottom red quantile removed (%) / 10')
    ylabel('Top Q quantile removed (%) / 10')
    title(sprintf('%s, %3.2f gof spread', time, range(tmp(:))))
end
exportgraphics(gcf, fullfile(save_direc, 'thermal_energies_gofs.png'), 'Resolution', 600)

figure('Position', [10, 60, 1500, 800])
tiledlayout(3,4)
for ie = 1:11
    nexttile
    hold on
    tmp = squeeze(thermal_energies_error(1:qlim,:,ie));
    imagesc(tmp)
    scatter(choice(2), choice(1), 20, 'r', 'filled')
    colorbar
    clim([0, 50])
    xlim([0.5,0.5+size(tmp,2)])
    ylim([0.5,0.5+size(tmp,1)])
    xlabel('Bottom red quantile removed (%) / 10')
    ylabel('Top Q quantile removed (%) / 10')
    title(sprintf('%s, %i eV spread', time, round(range(tmp(:)))))
end
exportgraphics(gcf, fullfile(save_direc, 'thermal_energies_error.png'), 'Resolution', 600)

figure('Position', [10, 60, 1500, 800])
tiledlayout(3,4)
for ie = 1:11
    nexttile
    hold on
    tmp = squeeze(thermal_energies_maxima(1:qlim,:,ie));
    imagesc(tmp)
    scatter(choice(2), choice(1), 10, 'r', 'filled')
    colorbar
    xlim([0.5,0.5+size(tmp,2)]); ylim([0.5,0.5+size(tmp,1)])
    xlabel('Bottom red quantile removed (%) / 10')
    ylabel('Top Q quantile removed (%) / 10')
    title(sprintf('%s, %i eV spread', time, round(range(tmp(:)))))
end
exportgraphics(gcf, fullfile(save_direc, 'thermal_energies_maxima.png'), 'Resolution', 600)

%%
for ie=1:11
    tmp = squeeze(thermal_energies_gofs(:,:,ie));
    [~,i] = max(tmp);
    [~,j] = max(max(tmp));
    fprintf('highest gof at (%i,%i) for event %i\n', i(j), j, ie)
end

for ie=1:11
    tmp = squeeze(thermal_energies_error(:,:,ie));
    [~,i] = min(tmp);
    [~,j] = min(min(tmp));
    fprintf('lowest error at (%i,%i) for event %i\n', i(j), j, ie)
end

%%
close all
for ie=1:11
    t = thermal_energies_fitted(choice(1),choice(2),ie);
    tmp = 100 * (squeeze(thermal_energies_fitted(:,:,ie)) / t - 1);
    histogram(tmp(:),linspace(-30,30,16))
    xlim([-30,30])
    input('aaaaa')
end