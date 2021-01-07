function res_plot(T, n, z, brd, tit, ticks, log_scale)
    figure;
    ax = axes;
    % dotsize = 100;  %adjust as needed
    x = n;
    y = T;
    
    z(z < brd(1)) = brd(1);
    z(z > brd(2)) = brd(2);
    % scatter3(x, y, z, dotsize, z, 'filled');
    surf(x, y, z);
    % title(tit, 'interpreter', 'none');
    title(tit, 'interpreter', 'tex');
    %title(tit);
    xlabel('n');
    ylabel('T');
    if(log_scale)
        set(ax,'XScale','log')
        set(ax,'YScale','log')
    end
    cb = colorbar;
    % ticks_labels = cellstr(num2str((log(ticks) ./ log(10))'));
    ticks_labels = cellstr(num2str(ticks'));
    set(cb, 'Ticks', log(ticks) ./ log(10));
    set(cb, 'TickLabels', ticks_labels);
    
    %view(0, -90);
    view(0, 90);
end

