function [names, colors] = def_names()
    names.frames_path = 'frames';
    names.E_filename = 'E.txt';
    names.T_filename = 'T.txt';
    names.time_filename = 'Time.txt';
    names.pressure_filename = 'Pressure.txt';
    names.head_filename = 'head.txt';
    names.params_filename = 'param.dat';
    names.data_path = 'DATA';
    names.frame_file_ext = 'xyz';
    names.x_clmn = 1;
    names.y_clmn = 2;
    names.z_clmn = 3;
    names.vx_clmn = 4;
    names.vy_clmn = 5;
    names.vz_clmn = 6;    
    
    colors = getMyColors;    
end

