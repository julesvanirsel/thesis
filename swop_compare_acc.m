direc = fullfile('data', 'paper2', 'dasc_data');
event_fn = fullfile('data', 'paper2', 'event_data.txt');
lines = readlines(event_fn)';

for li = 1:length(lines)-2
    l = lines(li+1);
    data = strsplit(l);
    time = datetime(data(2), 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
    te = str2double(data(13));
    time.Format = 'uuuuMMdd';

    close all
    figure('Position', [10, 50, 1500, 800])
    tiledlayout(2,4, 'TileSpacing', 'tight', 'Padding', 'tight')

    for is_acc = [1, 0]
        if is_acc
            h5fn = fullfile(direc, sprintf('INV_PKR_%s_%i_ACC.h5', time, second(time, 'secondofday')));
            row_title = sprintf('Accelerated Maxwellian (T_s = %i eV): ', te);
            energy_title = 'U_d (keV)';
        else
            h5fn = fullfile(direc, sprintf('INV_PKR_%s_%i.h5', time, second(time, 'secondofday')));
            row_title = 'Unaccelerated Maxwellian (T_s = E_0): ';
            energy_title = 'E_0 (keV)';
        end
        glat = h5read(h5fn, '/Coordinates/Latitude');
        glon = h5read(h5fn, '/Coordinates/Longitude');
        Q = h5read(h5fn, '/Derived/Energy/Flux');
        E0 = h5read(h5fn, '/Derived/Energy/Characteristic');
        SIGP = h5read(h5fn, '/Derived/Conductance/Pedersen');
        SIGH = h5read(h5fn, '/Derived/Conductance/Hall');
        % red = h5read(h5fn, '/Optical/Red');
        % blue = h5read(h5fn, '/Optical/Blue');
        % darkframes = h5read(h5fn, '/Optical/Preprocessing/Darkframes');
        % glat(isnan(E0)) = nan;
        % glon(isnan(E0)) = nan;
        % glat(red < quantile(red(:), 0.1)) = nan;
        E0(isnan(E0)) = 100; % otherwise will crash GEMINI
        E0(isnan(Q)) = nan;

        lim.x = [min(glon(not(isnan(Q)))), max(glon(not(isnan(Q))))];
        lim.y = [min(glat(not(isnan(Q)))), max(glat(not(isnan(Q))))];
        if is_acc
            lim.Q = [quantile(Q(:), 0.01), quantile(Q(:), 0.99)];
            lim.E0 = [quantile(E0(:), 0.01), quantile(E0(:), 0.99)];
            lim.SIG = [quantile([SIGP(:); SIGH(:)], 0.01), quantile([SIGP(:); SIGH(:)], 0.99)];
        end

        nexttile(is_acc*4+1)
        pcolor(glon, glat, Q)
        shading flat
        xlim(lim.x); ylim(lim.y); clim(lim.Q)
        xlabel('Geodetic Longitude (°)')
        ylabel('Geodetic Latitude (°)')
        colorbar
        title(sprintf('%s\nQ (mW/m^2)', row_title))

        nexttile(is_acc*4+2)
        pcolor(glon, glat, E0/1e3)
        shading flat
        xlim(lim.x); ylim(lim.y); clim(lim.E0/1e3)
        xlabel('Geodetic Longitude (°)')
        ylabel('Geodetic Latitude (°)')
        colorbar
        title(energy_title)

        nexttile(is_acc*4+3)
        pcolor(glon, glat, SIGP)
        shading flat
        xlim(lim.x); ylim(lim.y); clim(lim.SIG)
        xlabel('Geodetic Longitude (°)')
        ylabel('Geodetic Latitude (°)')
        colorbar
        title('\Sigma_P (S)')

        nexttile(is_acc*4+4)
        pcolor(glon, glat, SIGH)
        shading flat
        xlim(lim.x); ylim(lim.y); clim(lim.SIG)
        xlabel('Geodetic Longitude (°)')
        ylabel('Geodetic Latitude (°)')
        colorbar
        title('\Sigma_H (S)')
    end
    
    % input('Continue ... ')
    time.Format = 'uuuuMMdd';
    filename = sprintf('acc_comp_%s_%i.png', time, second(time, 'secondofday'));
    filename = fullfile('data', 'paper2', 'dasc_data', filename);
    exportgraphics(gcf, filename, 'Resolution', 600)

    direc_table = fullfile(direc, 'glow', sprintf('23%03i_%i', ...
        day(time, 'dayofyear'), second(time, 'secondofday')));
    table = jules.glow.read(direc_table);
    table_acc = jules.glow.read([direc_table, '_acc']);
    [~, iq] = min(abs(table.qvec - 10));
    evs = (1:3:10)*1e3;
    iecs = nan(size(evs));
    for i = 1:length(iecs)
        [~, iecs(i)] = min(abs(table.ecvec - evs(i)));
    end
    ialt = 56;
    alt = table.altvec(1:ialt);

    figure('Position', [10, 50, 600, 600])
    hold on
    
    for iec = iecs
        ne = squeeze(log10(table.ne(iq, iec, 1:ialt)));
        ne_acc = squeeze(log10(table_acc.ne(iq, iec, 1:ialt)));
        clr = [(iec-iecs(1))/range(iecs), 0, (iecs(end)-iec)/range(iecs)];
        qval = table.qvec(iq);
        ecval = table.ecvec(iec);
        plot(ne, alt, 'Color', clr, 'DisplayName', ...
            sprintf('Q = %.0f mW/m^2, E_0 = %.0f keV', qval, ecval/1e3))
        plot(ne_acc, alt, '--', 'Color', clr, 'DisplayName', ...
            sprintf('Q = %.0f mW/m^2, T_s = %.0f keV', qval, ecval/1e3))
    end
    xlabel('log_{10} Electron density (m^{-3})')
    ylabel('Altitude (km)')
    legend('Location', 'northwest')

    filename = sprintf('acc_comp_ne_%s_%i.png', time, second(time, 'secondofday'));
    filename = fullfile('data', 'paper2', 'dasc_data', filename);
    exportgraphics(gcf, filename, 'Resolution', 600)
end
