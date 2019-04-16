function leg = set_plot(ax, Xlbl, Ylbl, tit, x_scale, y_scale)
    set(ax, 'XScale', x_scale);
    set(ax, 'YScale', y_scale);
    grid(ax, 'on');
    xlabel(ax, Xlbl);
    ylabel(ax, Ylbl);
    title(ax, tit);
    leg = legend(ax, 'show');
    %set(leg, 'Interpreter', 'none');
    set(leg, 'Location', 'best');    
end

