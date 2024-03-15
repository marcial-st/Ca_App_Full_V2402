function [NCLS] = classGetDistance(NCLS,CELLS,inj_band)

cells_online = CELLS(contains({CELLS.status},'ONLINE'));            % Get only ONLINE cells

class_number = max(size(NCLS));

f = waitbar(0,'Thanks for waiting, distance computation in progress...');

for class_i = 1:1:class_number    
    cells_class_i = cells_online(NCLS(class_i).cVector);      % Filter only Class 2 cells    
    
    n_prof = max(size(cells_class_i));
    cells_distances = zeros(1,n_prof);
    for i=1:1:n_prof
        if (round(cells_class_i(i).xy_est(1)) >= inj_band(end))
            cells_class_i(i).inj_dist = round(cells_class_i(i).xy_est(1))-inj_band(end);
        elseif (round(cells_class_i(i).xy_est(1)) <= inj_band(1))
            cells_class_i(i).inj_dist = inj_band(1)-round(cells_class_i(i).xy_est(1));
        else
            cells_class_i(i).inj_dist = 0;
        end
        cells_distances(i) = cells_class_i(i).inj_dist;
    end    
    
    if mean(cells_distances)==0
        NCLS(class_i).class_dist_mean = 1;
    else
        NCLS(class_i).class_dist_mean = mean(cells_distances);
    end
    
    if std(cells_distances)==0
        NCLS(class_i).class_dist_std = 1;
    else
        NCLS(class_i).class_dist_std = std(cells_distances);
    end        
    NCLS(class_i).class_dist_max = max(cells_distances);
    waitbar(class_i/class_number,f);
end
close(f)