function [global_cells,cells_offline2map] = get_offline_cells(delta_frames,exp_data,Ts)
en_add_offline_cells = 1;
delta_frames = 0;

CELLS = exp_data.CELLS;
cells_online = CELLS(contains({CELLS.status},'ONLINE'));            % Get only ONLINE cells
cells_offline = CELLS(contains({CELLS.status},'OFFLINE'));

if en_add_offline_cells
    i_2map = 1;
    f_min = floor(delta_frames+exp_data.params.tstim/Ts);
    n_frames = length(cells_online(1).ca_profile);
    disp(strcat("looking for offline profiles beyond ",num2str(f_min)," frames"))
    for i = 1:1:length(cells_offline)
        if length(cells_offline(i).ca_profile)>f_min
            padval = cells_offline(i).ca_profile(end);
            i_frame = length(cells_offline(i).ca_profile);
            cells_offline(i).ca_profile(i_frame+1:n_frames) = padval;
            cells_offline2map(i_2map) = cells_offline(i); 
            i_2map = i_2map+1;
        end
    end
    global_cells = [cells_online,cells_offline2map];    
end      