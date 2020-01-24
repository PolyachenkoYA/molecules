function draw_diff(D_th, D_int, D_Einst, time, r2, lin_fit, v_cor, cubo_time_cut, scale, to_draw)
    if(~exist('scale','var'))
        scale = 'linear';
    end
    if(~exist('to_draw','var'))
        to_draw = [1 1 1];
    end    
    if(length(to_draw) < 3)
        to_draw = [to_draw ones(1, 3 - length(to_draw))];
    end
    cubo_i = time < cubo_time_cut;
    time = time(cubo_i);
    r2 = r2(cubo_i);
    v_cor = v_cor(cubo_i);
    
    if(to_draw(1))
        [fig, ax, leg] = getFig('time', '$\langle r^2 \rangle$',...
            ['$D_{Einst} = ' num2str(D_Einst) ', D_{th} = ' num2str(D_th) ' | err = ' num2str(eps_err(D_Einst, D_th) * 100) ' \%$'],...
            scale, scale);
        plot(time, r2, 'DisplayName', 'exp');
        plot(time, polyval(lin_fit, time), '--', 'DisplayName', 'line fit');
        %plot(time, r2 ./ (time * 6 * D_th), 'DisplayName', 'exp');
    end

    if(to_draw(2))
        [fig2, ax2, leg2] = getFig('time', '$\langle v(0)v(t) \rangle$',...
            ['$D_{cubo} = ' num2str(D_int) ', D_{th} = ' num2str(D_th) ' | err = ' num2str(eps_err(D_int, D_th) * 100) ' \%$'],...
                                   'linear', 'linear');
        plot(time, v_cor, 'DisplayName', 'exp');
        %plot(time, polyval(lin_fit, time), '--', 'DisplayName', 'line fit');
    end

    if(to_draw(3))
        [fig3, ax3, leg3] = getFig('time', '|\langle v(0)v(t) \rangle|',...
            ['$D_{cubo} = ' num2str(D_int) ', D_{th} = ' num2str(D_th) ' | err = ' num2str(eps_err(D_int, D_th) * 100) ' \%$'],...
                                   'linear', 'log');
        plot(time, abs(v_cor), 'DisplayName', 'exp');
    end
end

