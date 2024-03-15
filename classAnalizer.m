function NCLS = classAnalizer(NCLS,CELLS,en_fit_analysis)

% cells_online = CELLS(contains({CELLS.status},'ONLINE'));            % Get only ONLINE cells
cells_online = CELLS;
n_online = max(size(cells_online));

class_number = max(size(NCLS));

f = waitbar(0,'Please wait...');

for class_i = 1:1:class_number    
    cells_class_i = cells_online(NCLS(class_i).cVector);      % Filter only Class 2 cells    
    n_prof_i = max(size(cells_class_i));
    n_frames_i = max(size(cells_class_i(1).ca_profile));            % only valid since Class2 are online cells
    t = 1:1:n_frames_i;

    profiles_i = reshape([cells_class_i(:).ca_profile],n_frames_i,[]); % profiles array Profile-Column
    
    NCLS(class_i).cPercentage = n_prof_i/n_online;
    
    if en_fit_analysis
        NCLS(class_i).fitParams = classFitStats(profiles_i,20);
        NCLS(class_i).prof_fit  = NCLS(class_i).fitParams(1,1).*exp(-((t-NCLS(class_i).fitParams(2,1))/NCLS(class_i).fitParams(3,1)).^2) + NCLS(class_i).fitParams(4,1).*exp(-((t-NCLS(class_i).fitParams(5,1))/NCLS(class_i).fitParams(6,1)).^2) + NCLS(class_i).fitParams(7,1).*exp(-((t-NCLS(class_i).fitParams(8,1))/NCLS(class_i).fitParams(9,1)).^2);
    end

    bar_message = strcat("Extracting stats per class ",num2str(class_i));
    waitbar(class_i/class_number,f,bar_message);
end

close(f)