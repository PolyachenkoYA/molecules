function res_plot(ax, T, n, z, brd, tit, ticks, log_scale)
    % x=n, y=T & z' because axis look wrong on the final picture
    x = n;
    y = T;
    
    z(z < brd(1)) = brd(1);
    z(z > brd(2)) = brd(2);
    surf(ax, x, y, z, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    title(ax, tit, 'interpreter', 'tex');
    xlabel(ax, 'n (L-J)');
    ylabel(ax, 'T (L-J)');
    a = 0.01;
    xlim(ax, [(min(n) - max(n)*a) max(n)]);
    ylim(ax, [(min(T) - max(T)*a) max(T)]);
    if(log_scale)
        set(ax,'XScale','log')
        set(ax,'YScale','log')
    end
    cb = colorbar(ax);
    %caxis(brd);
    % ticks_labels = cellstr(num2str((log(ticks) ./ log(10))'));
    ticks_labels = cellstr(num2str(ticks'));
    set(cb, 'Ticks', log(ticks) ./ log(10));
    set(cb, 'TickLabels', ticks_labels);
    
    %view(ax, 0, -90);
    view(ax, 0, 90);
end

