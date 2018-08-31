function plot_lmd(T,n,z, brd, tit, log_scale)
    figure;
    ax = axes;
    dotsize = 100;  %adjust as needed
    x = n;
    y = T;    
    
    z(z < brd(1)) = brd(1);
    z(z > brd(2)) = brd(2);
    scatter3(x, y, z, dotsize, z, 'filled');
%     if(isempty(brd))
%         z(z < -0.5) = -0.5;
%         z(z > 0.5) = 0.5;        
%         scatter3(x, y, z, dotsize, z, 'filled');
%     else
%         scatter3(x, y, z, dotsize, z > brd, 'filled');
%         tit = [tit ' | dark == err < ' num2str(10^brd)];
%     end
    title(tit, 'interpreter', 'none');
    xlabel('n');
    ylabel('T');
    if(log_scale)
        set(ax,'XScale','log')
        set(ax,'YScale','log')
    end
    colorbar;
    %view(0, -90);
    view(0, 90);
end

