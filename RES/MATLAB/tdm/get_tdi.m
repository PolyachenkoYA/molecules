function td_i = get_tdi(v2_err, a)
    td_i = length(v2_err);
    while(v2_err(td_i) < a)
        td_i = td_i - 1;        
        if(td_i <= 0)
            disp('ERROR: td not found');
            td_i = 1;
            return;
        end
    end
end

