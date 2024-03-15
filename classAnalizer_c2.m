function [NCLS] = classAnalizer_c2(NCLS,cells_online,ca_prof_smooth,params,inj_band,pix2um,global_events,global_events_bin,fitmodel_c2_basal,fitmodel_c2_response,th_vals,en_ignore_low_resp)

n_online = max(size(cells_online));
Ts = 3; % TODO: Connect Ts to highest level in case user change this value from GUI
window_f = 30;
respath = "C:\Users\marci\OneDrive\Documentos\MATLAB\Results";
if params.tstart==0
    frame_start = 1;    
else
    frame_start = floor(params.tstart/Ts);
end
frame_inj = floor(params.tstim/Ts);
frame_end = floor(params.tend/Ts);


band_offset =  floor((inj_band(end)-inj_band(1))/2);
en_fit_analysis = 1;
en_plot_event = 0;
fitmodel = fitmodel_c2_response;

class_number = max(size(NCLS));

fwb = waitbar(0,'Please wait...');

if class_number==4
    disp("TODO: code or call thresholding")
elseif class_number==6
    disp("TODO: code or call thresholding")
elseif class_number==8
    disp("TODO: ok continue, nothing else to do")
end

for class_i = 3
    idx_vector = NCLS(class_i).cVector;
    cells_class_i = cells_online(NCLS(class_i).cVector);      % Filter only Class 2 cells    
    n_prof_i = max(size(cells_class_i));
    n_cells_i = length(cells_class_i);
    n_frames_i = max(size(cells_class_i(1).ca_profile));            % only valid since Class2 are online cells
    t = (0:1:n_frames_i-1)*Ts;

    profiles_i = reshape([cells_class_i(:).ca_profile],n_frames_i,[])'; 
    global_events_i     = [global_events(NCLS(class_i).cVector,:)];
    global_events_bin_i = [global_events_bin(NCLS(class_i).cVector,:)];
     
    basal_lines = mean(profiles_i(:,frame_start:frame_inj),2);
    basal_level.mean = mean(basal_lines);
    basal_level.std  = std(basal_lines);
    % 

    %% Plateau & Ratio stats
    th_plateu = 0.2;
    delta_plateau_ratio = (1-th_plateu)/3;
    [cplateau_1,cplateau_0,c_vector_plateau,plateau_ratio] = cmodPlateu2(ca_prof_smooth(idx_vector,:),frame_inj,frame_start,frame_inj,th_plateu,1,window_f);
    plateau_ratio    
    prob_plateu = sum(c_vector_plateau)/length(c_vector_plateau)
    NCLS(class_i).prob_plateau = prob_plateu;
    
    plateau_range_low    = th_plateu+delta_plateau_ratio;
    plateau_range_medium = th_plateu+2*delta_plateau_ratio;
    plateau_range_high   = 1;
    
    plateau_prob_low    = sum(and(plateau_ratio>=th_plateu,plateau_ratio<plateau_range_low))/sum(c_vector_plateau);
    plateau_prob_medium = sum(and(plateau_ratio>=plateau_range_low,plateau_ratio<plateau_range_medium))/sum(c_vector_plateau);
    plateau_prob_high   = sum(and(plateau_ratio>=plateau_range_medium,plateau_ratio<=plateau_range_high))/sum(c_vector_plateau);  
    if (plateau_prob_low + plateau_prob_medium + plateau_prob_high) ~= 1  % plateau_prob_check
        error('classAnalizer_c2 (67): plateau_prob_check failed')
    end
    
    NCLS(class_i).th_plateu = th_plateu;        
    NCLS(class_i).plateau_range_low    = plateau_range_low;    
    NCLS(class_i).plateau_range_medium = plateau_range_medium;
    NCLS(class_i).plateau_range_high   = plateau_range_high;
    
    NCLS(class_i).plateau_prob_low    = plateau_prob_low;    
    NCLS(class_i).plateau_prob_medium = plateau_prob_medium;
    NCLS(class_i).plateau_prob_high   = plateau_prob_high;
    
    c_vector_secpeak = sum(global_events_bin_i(:,frame_inj+window_f:frame_end),2)>0;
    secpeak_npeak = sum(global_events_bin_i(:,frame_inj+window_f:frame_end),2)
    prob_secpeak = sum(c_vector_secpeak)/length(c_vector_secpeak)
    NCLS(class_i).prob_secpeak = prob_secpeak;
    
    c_vector_palteau_secpeak = bitand(c_vector_plateau,c_vector_secpeak);
    prob_plateau_secpeak = sum(c_vector_palteau_secpeak)/length(c_vector_palteau_secpeak)
    NCLS(class_i).prob_plateau_secpeak = prob_plateau_secpeak;

   %% Computation: 
   % basal level per distance
   % first peak prominence

    i_cell_r = 0; i_cell_l = 0; cells_inside_inj=0;
    cells_ci_loc_x_l = [];cells_ci_loc_x_r = []; 
    cells_ci_d2inj_l_b = []; cells_ci_d2inj_r_b = []; cells_ci_d2inj = [];  % X
    cells_ci_loc_y_l = []; cells_ci_loc_y_r = []; % Y 
    basal_level_r = []; basal_level_l = [];
    main_peak_prom = []; main_peak_prom_r = [];   main_peak_prom_l = [];
    main_peak_time = []; main_peak_time_r = [];   main_peak_time_l = [];
    idx_rigth = []; idx_left = [];
    n_cells_i
    for i_cell = 1:1:n_cells_i
            cells_ci_loc(i_cell,1:2) = cells_class_i(i_cell).xy_hist(1,:);
            if not(isempty(inj_band))
                locs_all = [];
                peaks = find(global_events_bin_i(i_cell,frame_inj:frame_end)); % Get peaks from Event Matrix only after injury region              
                if ~isempty(peaks)
                    n_pks = length(peaks);        
                    peaks = peaks - 1;

                    if peaks(1) <= window_f  % Check if there is an event in the injury time window, 30 samples, 90s is the default vaue, it is hardcoded
                        [~,locs_all] = findpeaks(ca_prof_smooth(NCLS(class_i).cVector(i_cell),:)); % get all the signal peaks to fix the shifted ones due to 80/20 criteria         
                        locs_all = locs_all(locs_all>frame_inj+peaks(1));
                        if ((n_pks>1)&&(global_events(idx_vector(i_cell),frame_inj+peaks(1))<global_events(idx_vector(i_cell),frame_inj+peaks(2)))&&(peaks(2)<=(window_f))&&en_ignore_low_resp) % New criteria to filter artifacts right after the stimulus, the first peak is ignored if the second one has a greater prominence
                            main_peak_time(i_cell,:) = (locs_all(2)-frame_inj)*Ts;
                            main_peak_prom(i_cell,:) = global_events_i(i_cell,frame_inj+peaks(2));
                        else  % if en_ignore_response==0 then no event is shifted or ignored
                            main_peak_time(i_cell,:) = (locs_all(1)-frame_inj)*Ts;
                            main_peak_prom(i_cell,:) = global_events_i(i_cell,frame_inj+peaks(1));
                        end
                    else
                        disp(strcat("WARNING: PEAK OUT OF INJURY TIME WINDOW!!! ",num2str(i_cell)))
                        main_peak_time(i_cell,:) = peaks(1)*Ts
                        main_peak_prom(i_cell,:) = global_events_i(i_cell,frame_inj+peaks(1));
                        if global_events_i(i_cell,frame_inj+peaks(1))==0   % Check to trow a valid prominence value
                           error("Got an event with prominence value = 0")
                        end
                    end
                    if en_plot_event  % PLOT selected event on smoothed profile to visually validate by health exeperts
                        fig_c0 = figure;
                        plot(t,ca_prof_smooth(NCLS(class_i).cVector(i_cell),:));
                        hold on
                        plot(t(frame_inj+peaks(1)),ca_prof_smooth(NCLS(class_i).cVector(i_cell),frame_inj+peaks(1)),'r*');
                        if ~isempty(locs_all)
                           plot(t(locs_all),ca_prof_smooth(NCLS(class_i).cVector(i_cell),locs_all),'bv');
                        end
                        plot(t(frame_inj+floor(main_peak_time(i_cell,:)/Ts)),ca_prof_smooth(NCLS(class_i).cVector(i_cell),(frame_inj+floor(main_peak_time(i_cell,:)/Ts))),'m+','MarkerSize',10)
                        title(strcat("Class2 profile ",num2str(NCLS(class_i).cVector(i_cell))," i_cell=",num2str(i_cell), " event_index=",num2str(peaks(1))),'Interpreter','none')
                        hold off                       
                        ylim([0.5 1.5])
                        xline(frame_inj*Ts,'--')
                        yline(basal_lines(i_cell),'--')
                        saveas(fig_c0,fullfile(respath,strcat("IMG_B_",num2str(i_cell))),'png')                        
                        pause(0.01)
                    end
                end                   
                                                            %% Feature Extraction based on position
                                                            %
                if inj_band(1)>=cells_ci_loc(i_cell,1)      % Left side feature extraction
                    i_cell_l = i_cell_l+1;

                    cells_ci_loc_x_l(i_cell_l,:) = cells_ci_loc(i_cell,:); 
                    cells_ci_loc_y_l(i_cell_l,:) = cells_ci_loc(i_cell,2); 
                    cells_ci_d2inj(i_cell) = abs(inj_band(1)+band_offset-cells_ci_loc(i_cell,1))*pix2um;
                    basal_level_l(i_cell_l) =  mean(profiles_i(i_cell,frame_start:frame_inj));  
                    idx_left(i_cell_l) = i_cell;
                    if ~isempty(peaks)                        
                        main_peak_time_l(i_cell_l,:) = peaks(1)*Ts;
                        main_peak_prom_l(i_cell_l,:) = global_events_i(i_cell,frame_inj+peaks(1));
                    end                                      %%
                elseif inj_band(end)<=cells_ci_loc(i_cell,1) % Right side feature extraction
                    i_cell_r = i_cell_r+1;
                    cells_ci_loc_x_r(i_cell_r,:) = cells_ci_loc(i_cell,:); 
                    cells_ci_loc_y_r(i_cell_r,:) = cells_ci_loc(i_cell,2); 
                    cells_ci_d2inj(i_cell) = abs(inj_band(end)-band_offset-cells_ci_loc(i_cell,1))*pix2um;
                    basal_level_r(i_cell_r) =  mean(profiles_i(i_cell,frame_start:frame_inj));                    
                    idx_rigth(i_cell_r) = i_cell;
                    if ~isempty(peaks)
                        main_peak_time_r(i_cell_r,:) = peaks(1)*Ts;
                        main_peak_prom_r(i_cell_r,:) = global_events_i(i_cell,frame_inj+peaks(1));
                    end                      
                else
                    cells_ci_loc(i_cell,:) = [];
                    cells_ci_d2inj(i_cell) = 9999;
                    cells_inside_inj = cells_inside_inj + 1;
                    warning("Got a cell inside the injury")
                end
            end            
            waitbar(i_cell/n_cells_i,fwb,'Extracting stats Class 2');
    end
    close(fwb)

    cells_ci_d2inj_l_c = abs((inj_band(1)+band_offset)-cells_ci_loc_x_l(:,1))*pix2um;    % Left side distance to inj center
    cells_ci_d2inj_r_c = abs((inj_band(end)-band_offset)-cells_ci_loc_x_r(:,1))*pix2um; % Rigth side distance to inj center

    cells_ci_d2inj_l_b = abs((inj_band(1))-cells_ci_loc_x_l(:,1))*pix2um;    % Left side distance to inj boundary
    cells_ci_d2inj_r_b = abs(cells_ci_loc_x_r(:,1)-(inj_band(end)))*pix2um; % Rigth side distance to inj boundary    

    %% DATA CHECK STAGE
    if (length(main_peak_time_r)+length(main_peak_time_l))==(n_cells_i-cells_inside_inj)
        disp("Data Check 1/5 Passed: length main_peak_time_* and  n_cells_i")
    else 
        disp("(length(main_peak_time_r)+length(main_peak_time_l))==n_cells_i")
        length(main_peak_time_r)
        length(main_peak_time_l)
        n_cells_i
        cells_inside_inj
        error("Data Check 1/5 Failed: Lenght missmatch main_peak_time_* and n_cells_i")
    end

    if (length(main_peak_prom_r)+length(main_peak_prom_l))==(n_cells_i-cells_inside_inj)
        disp("Data Check 2/5 Passed: length main_peak_prom_* and  n_cells_i")
    else 
        disp("(length(main_peak_prom_r)+length(main_peak_prom_l))==n_cells_i")
        length(main_peak_prom_r)
        length(main_peak_prom_l)
        n_cells_i
        cells_inside_inj
        error("Data Check 2/5 Failed: Lenght missmatch main_peak_prom_* and n_cells_i")
    end

    if (length(basal_level_r)+length(basal_level_l))==(n_cells_i-cells_inside_inj)
        disp("Data Check 3/5 Passed: length basal_level_* and  n_cells_i")
    else
        disp("(length(basal_level_r)+length(basal_level_l))==n_cells_i")
        length(basal_level_r)
        length(basal_level_l)
        n_cells_i
        cells_inside_inj
        error("Data Check 3/5 Failed: Lenght missmatch basal_level_* and n_cells_i")        
    end

    if ((length(basal_level_r)==length(main_peak_time_r))&&((length(basal_level_r)==length(main_peak_prom_r))))
        disp("Data Check 4/5 Passed: right side number of features")
    else
        disp("((length(basal_level_r)==length(main_peak_time_r))&&((length(basal_level_r)==length(main_peak_prom_r))))")
        length(basal_level_r)
        length(main_peak_time_r)
        length(basal_level_r)
        length(main_peak_prom_r)
        error("Data Check 4/5 Failed: Lenght missmatch rigth side number of features")
    end

    if ((length(basal_level_l)==length(main_peak_time_l))&&((length(basal_level_l)==length(main_peak_prom_l))))
        disp("Data Check 5/5 Passed: left side number of features")
    else
        disp("((length(basal_level_l)==length(main_peak_time_l))&&((length(basal_level_l)==length(main_peak_prom_l))))")
        length(basal_level_l)
        length(main_peak_time_l)
        length(basal_level_l)
        length(main_peak_prom_l)
        error("Data Check 5/5 Failed: Lenght missmatch left side number of features")
    end
    %

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xmax = max([max(cells_ci_d2inj_l_b),max(cells_ci_d2inj_r_b)])+10;
    dist_nsample = 100;
    distv = linspace(0,xmax,dist_nsample);  
    model_list = {'poly1','poly2','poly3','exp1','exp2'};
    n_models = length(model_list);
    colorsm = rand(n_models,3);
    en_several_models = 0;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                Plots & Fitting: Basal Level
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitmodel_basal = 'poly1';    
    fig_a1 = figure;
    subplot(1,3,1)
        hold on
        scatter(cells_ci_d2inj_l_b,basal_level_l,'g')
        rsquare = [];                   % 2 report A
        rms = [];                       % 2 report A
        cl_ft = [];                     % 2 report A
        class_feature = 'C2_Basal_L';   % 2 report A        
        if en_several_models == 1                         
            for i_model = 1:1:n_models
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_l_b,basal_level_l',fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report A
                rms(i_model,1) = gof.rmse;          % 2 report A
                cl_ft = [cl_ft;class_feature];      % 2 report A                
                coeffs1 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs1,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report A
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report A            
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            legend(model_list)

            cfit = fit(cells_ci_d2inj_l_b,basal_level_l',fitmodel_basal);
            coeffs1 = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs1,distv,fitmodel_basal,0,0);
            plot(distv,fit_line_l,'m','LineWidth',1.5);
        else
            cfit = fit(cells_ci_d2inj_l_b,basal_level_l',fitmodel_basal);
            coeffs1 = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs1,distv,fitmodel_basal,0,0);
            plot(distv,fit_line_l,'m','LineWidth',1.5);
        end
        yline(0,'r--')
        title("Left Boundary")
        xlabel("Distance um")
        ylabel("Basal level")
        xlim([0 xmax])
        ylim([min(basal_lines)-0.1 max(basal_lines)+0.1])
        hold off
    subplot(1,3,2)
        hold on
        scatter(-cells_ci_d2inj_l_b,basal_level_l,'g')
        scatter(cells_ci_d2inj_r_b,basal_level_r,'g')
        % scatter(cells_ci_d2inj,basal_lines,'g')
        % cfit = fit(cells_ci_d2inj',basal_lines,'poly1');    
        % coeffs = coeffvalues(cfit);
        % fit_line = coeffs(1)*cells_ci_d2inj+coeffs(2);
        % NCLS(class_i).fitbasal_params = coeffs;
        % plot(cells_ci_d2inj,fit_line,'m',LineWidth=1.5);        
        title("Center")
        xlabel("Distance um")
        ylabel("Basal level")        
        xl = max([max(cells_ci_d2inj_l_b),max(cells_ci_d2inj_r_b)]);        
        xlim([-xl-10 xl+10]) 
        xline(0,'--')
        ylim([min(basal_lines)-0.1 max(basal_lines)+0.1])
        yline(0,'r--')
        hold off
    subplot(1,3,3)
        hold on
        scatter(cells_ci_d2inj_r_b,basal_level_r,'g')
        rsquare = [];                   % 2 report A
        rms = [];                       % 2 report A
        cl_ft = [];                     % 2 report A
        class_feature = 'C2_Basal_R';   % 2 report A
        model_list(1)=[];
         if en_several_models == 1 
            for i_model = 1:1:n_models
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_r_b,basal_level_r',fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report A
                rms(i_model,1) = gof.rmse;          % 2 report A
                cl_ft = [cl_ft;class_feature];      % 2 report A                
                coeffs2 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs2,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report A
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report A            
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            model_list{end+1} = 'Mean Fit';
            legend(model_list)

            cfit = fit(cells_ci_d2inj_r_b,basal_level_r',fitmodel_basal);
            coeffs2 = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs2,distv,fitmodel_basal,0,0);
            plot(distv,fit_line_r,'m','LineWidth',1.5);
        else
            cfit = fit(cells_ci_d2inj_r_b,basal_level_r',fitmodel_basal);
            coeffs2 = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs2,distv,fitmodel_basal,0,0);
            plot(distv,fit_line_r,'m','LineWidth',1.5);
        end            
        coeffs_m = mean([coeffs1;coeffs2],1);   % 2 model
        fit_line_m_basal = genProfileFit(coeffs_m,distv,fitmodel_basal,0,0);
        NCLS(class_i).fitparams_basal = coeffs_m;
        NCLS(class_i).fitparams_basal_model = fitmodel_basal;        
        plot(distv,fit_line_m_basal,'c','LineWidth',1.5)
        yline(0,'r--')
        title("Right Boundary")
        xlabel("Distance um")
        ylabel("Basal level")
        xlim([0 max([max(cells_ci_d2inj_l_b),max(cells_ci_d2inj_r_b)])+10])
        ylim([min(basal_lines)-0.1 max(basal_lines)+0.1])
        hold off
    sgtitle("Basal level Vs Distance to injury")
    saveas(fig_a1,fullfile(respath,"IMG_A1_basal"),'png')
    
    FeatureName = "C2_BasalLevel";                                                                          % 2 report B
    tableData = table(FeatureName,{fitmodel_basal},coeffs_m,'VariableNames',{'Feature','Model','Coeffs'});  % 2 report B
    writetable(tableData,'fitModels.txt','Delimiter','\t','WriteMode','append');                            % 2 report B   

    %%%% Thesis plot begin
        figure; 
        hold on;
        color_r = [180,180,180]/255;
        color_l = [180,180,180]/255;    
        color_m = [15,15,15]/255;
        marker_l = "o";
        marker_r = "+";
        size_axis = 20;
        size_label = 24;
        size_legend = 14;
        scatter(cells_ci_d2inj_l_b,basal_level_l,50,'MarkerEdgeColor',color_l,'Marker',marker_l) % Left data
        plot(distv,fit_line_l,'--',LineWidth=1.5,Color=color_l);        
        scatter(cells_ci_d2inj_r_b,basal_level_r,55,'MarkerEdgeColor',color_r,'Marker',marker_r) % Rigth data
        plot(distv,fit_line_r,'--',LineWidth=1.5,Color=color_r);        
        plot(distv,fit_line_m_basal,'-',LineWidth=1.5,Color=color_m)
        xlabel("Distance to injury (\mum)",'Interpreter','tex','FontSize',size_label)
        ylabel("Basal level (F_{340}/F_{380})",'Interpreter','tex','FontSize',size_label)
        set(gca,'FontSize',size_axis)
        xlim([0 max([max(cells_ci_d2inj_l_b),max(cells_ci_d2inj_r_b)])+10])
        % ylim([min([basal_level_l,basal_level_r])-0.1 max([basal_level_l,basal_level_r])+0.1])
        ylim([0.7 1.2])
        qw{1} = plot(nan,"o",'MarkerEdgeColor',color_l);
        qw{2} = plot(nan,"+",'MarkerEdgeColor',color_r);
        qw{3} = plot(nan, '-',LineWidth=1.5,Color=color_m);
        legend([qw{:}],{'Left side data','Rigth side data','Model'},'FontSize',size_legend)    
        hold off
    %%%% Thesis plot end 


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                Plots & Fitting: Peak Prominence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx_wf_l = find(main_peak_time_l<=window_f*Ts);
    idx_wf_r = find(main_peak_time_r<=window_f*Ts);

    clear model_list
    model_list = {'poly1','poly2','poly3','exp1','exp2'};
    fitmodel_prominence = 'exp1';
    fig_a1 = figure;    
    subplot(1,3,1)
        hold on
        scatter(cells_ci_d2inj_l_b(idx_wf_l),main_peak_prom_l(idx_wf_l),'g')
        rsquare = [];                   % 2 report A
        rms = [];                       % 2 report A
        cl_ft = [];                     % 2 report A
        class_feature = 'C2_PeakProm_L';% 2 report A        
        if en_several_models == 1                         
            for i_model = 1:1:n_models
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_l_b(idx_wf_l),main_peak_prom_l(idx_wf_l),fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report A
                rms(i_model,1) = gof.rmse;          % 2 report A
                cl_ft = [cl_ft;class_feature];      % 2 report A                
                coeffs1 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs1,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report A
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report A            
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            legend(model_list)
            cfit = []; coeffs1 = [];
            cfit = fit(cells_ci_d2inj_l_b(idx_wf_l),main_peak_prom_l(idx_wf_l),fitmodel_prominence);
            coeffs1 = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs1,distv,fitmodel_prominence,0,0);
            plot(distv,fit_line_l,'m','LineWidth',1.5);            
        else
            cfit = fit(cells_ci_d2inj_l_b(idx_wf_l),main_peak_prom_l(idx_wf_l),fitmodel_prominence);
            coeffs1 = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs1,distv,fitmodel_prominence,0,0);
            plot(distv,fit_line_l,'m','LineWidth',1.5);
        end
        yline(0,'r--')
        title("Left Boundary")
        xlabel("Distance um")
        ylabel("Peak Prominence")
        xlim([0 max([max(cells_ci_d2inj_l_c),max(cells_ci_d2inj_r_c)])+10])
        ylim([min(main_peak_prom)-0.1 max(main_peak_prom)+0.1])
        hold off
    subplot(1,3,2)
        hold on
        scatter(-cells_ci_d2inj_l_b(idx_wf_l),main_peak_prom_l(idx_wf_l),'g')
        scatter(cells_ci_d2inj_r_b(idx_wf_r),main_peak_prom_r(idx_wf_r),'g')
        % cfit = fit(cells_ci_d2inj',basal_lines,'poly1');    
        % coeffs = coeffvalues(cfit);
        % fit_line = coeffs(1)*cells_ci_d2inj+coeffs(2);
        % NCLS(class_i).fitbasal_params = coeffs;
        % plot(cells_ci_d2inj,fit_line,'m',LineWidth=1.5);
        yline(0,'r--')
        title("Center")
        xlabel("Distance um")
        ylabel("Peak Prominence")
        xl = max([max(cells_ci_d2inj_l_b(idx_wf_l)),max(cells_ci_d2inj_r_b(idx_wf_r))]);        
        xlim([-xl-10 xl+10]) 
        xline(0,'--')
        ylim([min(main_peak_prom)-0.1 max(main_peak_prom)+0.1])
        hold off
    subplot(1,3,3)
        hold on
        scatter(cells_ci_d2inj_r_b(idx_wf_r),main_peak_prom_r(idx_wf_r),'g')
        rsquare = [];                   % 2 report A
        rms = [];                       % 2 report A
        cl_ft = [];                     % 2 report A
        class_feature = 'C2_PeakProm_R';% 2 report A        
        model_list(1)=[];
         if en_several_models == 1 
            for i_model = 1:1:n_models
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_r_b(idx_wf_r),main_peak_prom_r(idx_wf_r),fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report A
                rms(i_model,1) = gof.rmse;          % 2 report A
                cl_ft = [cl_ft;class_feature];      % 2 report A                
                coeffs2 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs2,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report A
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report A            
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            model_list{end+1} = 'Mean Fit';
            legend(model_list)
            cfit = []; coeffs2 = [];
            cfit = fit(cells_ci_d2inj_r_b(idx_wf_r),main_peak_prom_r(idx_wf_r),fitmodel_prominence);
            coeffs2 = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs2,distv,fitmodel_prominence,0,0);
            plot(distv,fit_line_r,'m','LineWidth',1.5);            
        else
            cfit = fit(cells_ci_d2inj_r_b(idx_wf_r),main_peak_prom_r(idx_wf_r),fitmodel_prominence);
            coeffs2 = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs2,distv,fitmodel_prominence,0,0);
            plot(distv,fit_line_r,'m','LineWidth',1.5);
        end            
        coeffs_m = mean([coeffs1;coeffs2],1)   % 2 model
        fit_line_m_peakprom = genProfileFit(coeffs_m,distv,fitmodel_prominence,0,0);
        NCLS(class_i).fitparams_prominence = coeffs_m;
        NCLS(class_i).fitparams_prominence_model = fitmodel_prominence;        
        plot(distv,fit_line_m_peakprom,'c','LineWidth',1.5)
        yline(0,'r--')
        title("Right Boundary")
        xlabel("Distance um")
        ylabel("Peak Prominence")
        xlim([0 max([max(cells_ci_d2inj_l_c),max(cells_ci_d2inj_r_c)])+10])
        ylim([min(main_peak_prom)-0.1 max(main_peak_prom)+0.1])
        hold off
    sgtitle("Peak Prominence Vs Distance to injury")
    saveas(fig_a1,fullfile(respath,"IMG_A1_prominence"),'png')

    FeatureName = "C2_PeakProm";                                                                                                      % 2 report B
    tableData = table(FeatureName,{fitmodel_prominence},coeffs_m,-1/coeffs_m(2),'VariableNames',{'Feature','Model','Coeffs','Tau'});  % 2 report B
    writetable(tableData,'fitModels.txt','Delimiter','\t','WriteMode','append');                                                      % 2 report B   

    %%%% Thesis plot begin
        figure; 
        hold on;
        color_r = [180,180,180]/255;
        color_l = [180,180,180]/255;    
        color_m = [15,15,15]/255;
        marker_l = "o";
        marker_r = "+";
        size_axis = 20;
        size_label = 24;
        size_legend = 14;
        scatter(cells_ci_d2inj_l_b(idx_wf_l),main_peak_prom_l(idx_wf_l),50,'MarkerEdgeColor',color_l,'Marker',marker_l) % Left data
        plot(distv,fit_line_l,'--',LineWidth=1.5,Color=color_l);        
        scatter(cells_ci_d2inj_r_b(idx_wf_r),main_peak_prom_r(idx_wf_r),55,'MarkerEdgeColor',color_r,'Marker',marker_r) % Rigth data
        plot(distv,fit_line_r,'--',LineWidth=1.5,Color=color_r);        
        plot(distv,fit_line_m_peakprom,'-',LineWidth=1.5,Color=color_m)
        xlabel("Distance to injury (\mum)",'Interpreter','tex','FontSize',size_label)
        ylabel("MP (F_{340}/F_{380})",'Interpreter','tex','FontSize',size_label)
        set(gca,'FontSize',size_axis)
        xlim([0 max([max(cells_ci_d2inj_l_b(idx_wf_l)),max(cells_ci_d2inj_r_b(idx_wf_r))])+10])
        % ylim([min([main_peak_prom_l',main_peak_prom_r'])-0.1 max([main_peak_prom_l',main_peak_prom_r'])+0.1])
        ylim([-0.1 0.8])
        qw{1} = plot(nan,"o",'MarkerEdgeColor',color_l);
        qw{2} = plot(nan,"+",'MarkerEdgeColor',color_r);
        qw{3} = plot(nan, '-',LineWidth=1.5,Color=color_m);
        legend([qw{:}],{'Left side data','Rigth side data','Model'},'FontSize',size_legend)    
        hold off
    %%%% Thesis plot end 


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                Plots & Fitting: Peak Time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
    [~,idx_out_l] = rmoutliers(main_peak_time_l);
    [~,idx_out_r] = rmoutliers(main_peak_time_r);
    idx_out_l = not(idx_out_l);
    idx_out_r = not(idx_out_r);
    clear model_list
    model_list = {'poly1','poly2','poly3'};
    fitmodel_peaktime = 'poly1';    
    fig_a1 = figure;
    subplot(1,3,1)
        hold on
        scatter(cells_ci_d2inj_l_b(idx_wf_l),main_peak_time_l(idx_wf_l),'g')
        rsquare = [];                   % 2 report A
        rms = [];                       % 2 report A
        cl_ft = [];                     % 2 report A
        class_feature = 'C2_PeakTime_L';% 2 report A        
        if en_several_models == 1                         
            for i_model = 1:1:length(model_list)
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_l_b(idx_wf_l),main_peak_time_l(idx_wf_l),fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report A
                rms(i_model,1) = gof.rmse;          % 2 report A
                cl_ft = [cl_ft;class_feature];      % 2 report A                
                coeffs1 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs1,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end

            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            legend(model_list)

            cfit = fit(cells_ci_d2inj_l_b(idx_wf_l),main_peak_time_l(idx_wf_l),fitmodel_peaktime);
            coeffs1 = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs1,distv,fitmodel_peaktime,0,0);
            plot(distv,fit_line_l,'m','LineWidth',1.5);            
        else
            cfit = fit(cells_ci_d2inj_l_b(idx_wf_l),main_peak_time_l(idx_wf_l),fitmodel_peaktime);
            coeffs1 = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs1,distv,fitmodel_peaktime,0,0);
            plot(distv,fit_line_l,'m','LineWidth',1.5);
        end
        yline(0,'r--')
        title("Left Boundary")
        xlabel("Distance um")
        ylabel("Peak time")
        xlim([0 xmax])
        ylim([min(main_peak_time)-0.1 max(main_peak_time)+0.1])
        hold off
    subplot(1,3,2)
        hold on
        scatter(-cells_ci_d2inj_l_b(idx_wf_l),main_peak_time_l(idx_wf_l),'g')
        scatter(cells_ci_d2inj_r_b(idx_wf_r),main_peak_time_r(idx_wf_r),'g')
        % scatter(cells_ci_d2inj,main_peak_time,'g')
        % cfit = fit(cells_ci_d2inj',main_peak_time,'poly1');    
        % coeffs = coeffvalues(cfit);
        % fit_line = coeffs(1)*cells_ci_d2inj+coeffs(2);
        % NCLS(class_i).fitbasal_params = coeffs;
        % plot(cells_ci_d2inj,fit_line,'m',LineWidth=1.5);        
        yline(0,'r--')
        title("Center")
        xlabel("Distance um")
        ylabel("Peak time")        
        xl = max([max(cells_ci_d2inj_l_b(idx_wf_l)),max(cells_ci_d2inj_r_b(idx_wf_r))]);        
        xlim([-xl-10 xl+10]) 
        xline(0,'--')
        ylim([min(main_peak_time)-0.1 max(main_peak_time)+0.1])
        hold off
    subplot(1,3,3)
        hold on
        scatter(cells_ci_d2inj_r_b(idx_wf_r),main_peak_time_r(idx_wf_r),'g')
        rsquare = [];                   % 2 report A
        rms = [];                       % 2 report A
        cl_ft = [];                     % 2 report A
        class_feature = 'C2_PeakTime_R';% 2 report A        
        model_list(1)=[];
         if en_several_models == 1 
            for i_model = 1:1:length(model_list)
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_r_b(idx_wf_r),main_peak_time_r(idx_wf_r),fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report A
                rms(i_model,1) = gof.rmse;          % 2 report A
                cl_ft = [cl_ft;class_feature];      % 2 report A                
                coeffs2 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs2,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report A
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report A            
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            model_list{end+1} = 'Mean Fit';
            legend(model_list)

            cfit = fit(cells_ci_d2inj_r_b(idx_wf_r),main_peak_time_r(idx_wf_r),fitmodel_peaktime);
            coeffs2 = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs2,distv,fitmodel_peaktime,0,0);
            plot(distv,fit_line_r,'m','LineWidth',1.5);            
        else
            cfit = fit(cells_ci_d2inj_r_b(idx_wf_r),main_peak_time_r(idx_wf_r),fitmodel_peaktime);
            coeffs2 = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs2,distv,fitmodel_peaktime,0,0);
            plot(distv,fit_line_r,'m','LineWidth',1.5);
        end            
        coeffs_m = mean([coeffs1;coeffs2],1);   % 2 model
        fit_line_m_peaktime = genProfileFit(coeffs_m,distv,fitmodel_peaktime,0,0);
        NCLS(class_i).fitparams_peaktime = coeffs_m;
        NCLS(class_i).fitparams_peaktime_model = fitmodel_peaktime;        
        plot(distv,fit_line_m_peaktime,'c','LineWidth',1.5)
        yline(0,'r--')
        title("Right Boundary")
        xlabel("Distance um")
        ylabel("Peak time")
        xlim([0 max([max(cells_ci_d2inj_l_b(idx_wf_l)),max(cells_ci_d2inj_r_b(idx_wf_r))])+10])
        ylim([min(main_peak_time)-0.1 max(main_peak_time)+0.1])
        hold off
    sgtitle("Peak time Vs Distance to injury")
    saveas(fig_a1,fullfile(respath,"IMG_A1_peaktime"),'png')

    FeatureName = "C2_PeakTime";                                                                              % 2 report B
    tableData = table(FeatureName,{fitmodel_peaktime},coeffs_m,'VariableNames',{'Feature','Model','Coeffs'}); % 2 report B
    writetable(tableData,'fitModels.txt','Delimiter','\t','WriteMode','append');                              % 2 report B

    %%%% Thesis plot begin
        figure; 
        hold on;
        color_r = [180,180,180]/255;
        color_l = [180,180,180]/255;    
        color_m = [15,15,15]/255;
        marker_l = "o";
        marker_r = "+";
        size_axis = 20;
        size_label = 24;
        size_legend = 14;
        scatter(cells_ci_d2inj_l_b(idx_wf_l),main_peak_time_l(idx_wf_l),50,'MarkerEdgeColor',color_l,'Marker',marker_l) % Left data
        plot(distv,fit_line_l,'--',LineWidth=1.5,Color=color_l);        
        scatter(cells_ci_d2inj_r_b(idx_wf_r),main_peak_time_r(idx_wf_r),55,'MarkerEdgeColor',color_r,'Marker',marker_r) % Rigth data
        plot(distv,fit_line_r,'--',LineWidth=1.5,Color=color_r);        
        plot(distv,fit_line_m_peaktime,'-',LineWidth=1.5,Color=color_m)
        xlabel("Distance to injury (\mum)",'Interpreter','tex','FontSize',size_label)
        ylabel("MP Time (seconds)",'Interpreter','tex','FontSize',size_label)
        set(gca,'FontSize',size_axis)
        xlim([0 max([max(cells_ci_d2inj_l_b(idx_wf_l)),max(cells_ci_d2inj_r_b(idx_wf_r))])+10])
        ylim([0 100])
        % ylim([min([main_peak_time_l(idx_out_l)',main_peak_time_r(idx_out_r)'])-0.1 max([main_peak_time_l(idx_out_l)',main_peak_time_r(idx_out_r)'])+0.1])
        qw{1} = plot(nan,"o",'MarkerEdgeColor',color_l);
        qw{2} = plot(nan,"+",'MarkerEdgeColor',color_r);
        qw{3} = plot(nan, '-',LineWidth=1.5,Color=color_m);
        legend([qw{:}],{'Left side data','Rigth side data','Model'},'FontSize',size_legend)    
        hold off
    %%%% Thesis plot end     

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                Plots & Fitting: Plateau
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    clear model_list
    model_list = {'poly1','poly2'};    
    fitmodel_plateau = 'poly1';    
    fig_a1 = figure;
    subplot(1,3,1)
        hold on        
        idx_aux_l = find(plateau_ratio(idx_left)>0);
        scatter(cells_ci_d2inj_l_b(idx_aux_l),plateau_ratio(idx_left(idx_aux_l)),'g')
        rsquare = [];                   % 2 report A
        rms = [];                       % 2 report A
        cl_ft = [];                     % 2 report A
        class_feature = 'C2_Plateau_L';% 2 report A        
        if en_several_models == 1                         
            for i_model = 1:1:length(model_list)
                i_model
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_l_b(idx_aux_l),plateau_ratio(idx_left(idx_aux_l)),fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report A
                rms(i_model,1) = gof.rmse;          % 2 report A
                cl_ft = [cl_ft;class_feature];      % 2 report A                
                coeffs1 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs1,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report A
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report A            
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            legend(model_list)

            cfit = fit(cells_ci_d2inj_l_b(idx_aux_l),plateau_ratio(idx_left(idx_aux_l)),fitmodel_plateau);
            coeffs1 = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs1,distv,fitmodel_plateau,0,0);
            plot(distv,fit_line_l,'m','LineWidth',1.5);
        else
            cfit = fit(cells_ci_d2inj_l_b(idx_aux_l),plateau_ratio(idx_left(idx_aux_l)),fitmodel_plateau);
            coeffs1 = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs1,distv,fitmodel_plateau,0,0);
            plot(distv,fit_line_l,'m','LineWidth',1.5);
        end
        yline(0,'r--')
        title("Left Boundary")
        xlabel("Distance um")
        ylabel("Plateau Ratio")
        xlim([0 xmax])
        ylim([min(plateau_ratio)-0.1 max(plateau_ratio)+0.1])
        hold off
    subplot(1,3,2)
        hold on
        idx_aux_l = find(plateau_ratio(idx_left)>0);
        scatter(-cells_ci_d2inj_l_b(idx_aux_l),plateau_ratio(idx_left(idx_aux_l)),'g')
        idx_aux_r = find(plateau_ratio(idx_rigth)>0);
        scatter( cells_ci_d2inj_r_b(idx_aux_r),plateau_ratio(idx_rigth(idx_aux_r)),'g')
        % scatter(cells_ci_d2inj,basal_lines,'g')
        % cfit = fit(cells_ci_d2inj',basal_lines,'poly1');    
        % coeffs = coeffvalues(cfit);
        % fit_line = coeffs(1)*cells_ci_d2inj+coeffs(2);
        % NCLS(class_i).fitbasal_params = coeffs;
        % plot(cells_ci_d2inj,fit_line,'m',LineWidth=1.5);        
        title("Center")
        xlabel("Distance um")
        ylabel("Plateau Ratio")        
        xl = max([max(cells_ci_d2inj_l_b),max(cells_ci_d2inj_r_b)]);        
        xlim([-xl-10 xl+10]) 
        xline(0,'--')
        ylim([min(plateau_ratio)-0.1 max(plateau_ratio)+0.1])
        yline(0,'r--')
        hold off
    subplot(1,3,3)
        hold on
        idx_aux_r = find(plateau_ratio(idx_rigth)>0);
        scatter(cells_ci_d2inj_r_b(idx_aux_r),plateau_ratio(idx_rigth(idx_aux_r)),'g')
        rsquare = [];                   % 2 report A
        rms = [];                       % 2 report A
        cl_ft = [];                     % 2 report A
        class_feature = 'C2_Plateau_R'; % 2 report A        
        model_list(1)=[];
         if en_several_models == 1 
            for i_model = 1:1:length(model_list)
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_r_b(idx_aux_r),plateau_ratio(idx_rigth(idx_aux_r)),fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report A
                rms(i_model,1) = gof.rmse;          % 2 report A
                cl_ft = [cl_ft;class_feature];      % 2 report A                
                coeffs2 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs2,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report A
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report A            
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            model_list{end+1} = 'Mean Fit';
            legend(model_list)

            cfit = fit(cells_ci_d2inj_r_b(idx_aux_r),plateau_ratio(idx_rigth(idx_aux_r)),fitmodel_plateau);
            coeffs2 = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs2,distv,fitmodel_plateau,0,0);
            plot(distv,fit_line_r,'m','LineWidth',1.5);
        else
            cfit = fit(cells_ci_d2inj_r_b(idx_aux_r),plateau_ratio(idx_rigth(idx_aux_r)),fitmodel_plateau);
            coeffs2 = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs2,distv,fitmodel_plateau,0,0);
            plot(distv,fit_line_r,'m','LineWidth',1.5);
        end            
        coeffs_m = mean([coeffs1;coeffs2],1);   % 2 model
        fit_line_m_plateau = genProfileFit(coeffs_m,distv,fitmodel_plateau,0,0);
        NCLS(class_i).fitparams_plateau = coeffs_m;
        NCLS(class_i).fitparams_plateau_model = fitmodel_plateau;        
        plot(distv,fit_line_m_plateau,'c','LineWidth',1.5)
        yline(0,'r--')
        title("Right Boundary")
        xlabel("Distance um")
        ylabel("Plateau Ratio")
        xlim([0 max([max(cells_ci_d2inj_l_b),max(cells_ci_d2inj_r_b)])+10])
        ylim([min(plateau_ratio)-0.1 max(plateau_ratio)+0.1])
        hold off
    sgtitle("Plateau Ratio Vs Distance to injury")
    saveas(fig_a1,fullfile(respath,"IMG_A1_plateau"),'png')

    FeatureName = "C2_Plateau";                                                                             % 2 report B
    tableData = table(FeatureName,{fitmodel_plateau},coeffs_m,'VariableNames',{'Feature','Model','Coeffs'});% 2 report B
    writetable(tableData,'fitModels.txt','Delimiter','\t','WriteMode','append');                            % 2 report B
    
    %%%% Thesis plot begin
        thb = 0.2;
        delta=(1-thb)/3;
        m1 = 0.2; m2 = m1+delta; m3=m2+delta;

        figure; 
        hold on;
        color_r = [180,180,180]/255;
        color_l = [180,180,180]/255;    
        color_m = [15,15,15]/255;
        marker_l = "o";
        marker_r = "+";
        size_axis = 20;
        size_label = 24;
        size_legend = 14;
        scatter(cells_ci_d2inj_l_b(idx_aux_l),plateau_ratio(idx_aux_l),50,'MarkerEdgeColor',color_l,'Marker',marker_l) % Left data
        plot(distv,fit_line_l,'--',LineWidth=1.5,Color=color_l);        
        scatter(cells_ci_d2inj_r_b(idx_aux_r),plateau_ratio(idx_aux_r),55,'MarkerEdgeColor',color_r,'Marker',marker_r) % Rigth data
        plot(distv,fit_line_r,'--',LineWidth=1.5,Color=color_r);        
        plot(distv,fit_line_m_plateau,'-',LineWidth=1.5,Color=color_m)
        xlabel("Distance to injury (\mum)",'Interpreter','tex','FontSize',size_label)
        ylabel("Plateau (ratio)",'Interpreter','tex','FontSize',size_label)
        set(gca,'FontSize',size_axis)
        xlim([0 max([max(cells_ci_d2inj_l_b(idx_aux_l)),max(cells_ci_d2inj_r_b(idx_aux_r))])+10])
        % ylim([min([plateau_ratio(idx_aux_l)',plateau_ratio(idx_aux_r)'])-0.1 max([plateau_ratio(idx_aux_l)',plateau_ratio(idx_aux_r)'])+0.1])
        ylim([-0.1 1.1])
        yline(m1,'-.','th1','Fontsize',12,'LabelVerticalAlignment','bottom')
        yline(m2,'-.','th2','Fontsize',12,'LabelVerticalAlignment','bottom')
        yline(m3,'-.','th3','Fontsize',12,'LabelVerticalAlignment','bottom')
        qw{1} = plot(nan,"o",'MarkerEdgeColor',color_l);
        qw{2} = plot(nan,"+",'MarkerEdgeColor',color_r);
        qw{3} = plot(nan, '-',LineWidth=1.5,Color=color_m);
        legend([qw{:}],{'Left side data','Rigth side data','Model'},'FontSize',size_legend)    
        hold off
    %%%% Thesis plot end     

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                Plots & Fitting: Summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    fig_a1 = figure;
    hold on
    plot(distv,fit_line_m_basal,'LineWidth',1.5)
    plot(distv,fit_line_m_peakprom,'LineWidth',1.5)
    plot(distv,fit_line_m_peaktime,'LineWidth',1.5)
    plot(distv,fit_line_m_plateau,'LineWidth',1.5)
    %% Get Range for three features    
    lim_basal    = find(fit_line_m_basal(1:end-1)>0 & fit_line_m_basal(2:end) < 0);
    lim_peakprom = find(fit_line_m_peakprom(1:end-1)>0 & fit_line_m_peakprom(2:end) < 0);
    lim_peaktime = find(fit_line_m_peaktime(1:end-1)>0 & fit_line_m_peaktime(2:end) < 0);
    lim_plateau = find(fit_line_m_plateau(1:end-1)>0 & fit_line_m_plateau(2:end) < 0);
    lim_global = min([lim_basal,lim_peakprom,lim_peaktime,lim_plateau]);

    if ~isempty(lim_basal)
        xline(lim_basal,'g--')
    end
    if ~isempty(lim_peakprom)
        xline(lim_peakprom,'g--')
    end
    if ~isempty(lim_peaktime)
        xline(lim_peaktime,'g--')
    end
    if ~isempty(lim_plateau)
        xline(lim_plateau,'g--')
    end
    if ~isempty(lim_global)
        xline(lim_global,'r--')
    end

    yline(0,'r--')
    legend({'Basal','Peak Prominence','Peak Time','Plateau'})
    hold off
    sgtitle("Summary Features Vs Distance to injury")
    xlabel("Distance um")
    ylabel("Magnitude")
    saveas(fig_a1,fullfile(respath,"IMG_A1_summary"),'png')

    % genReport(respath,strcat("Exp. 100528_4"),[],[])


    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    figure
    subplot(1,2,1)
        hold on
        scatter(cells_ci_loc_y_l*pix2um,basal_level_l,'c')
        cfit = fit(cells_ci_loc_y_l*pix2um,basal_level_l','poly1');    
        coeffs = coeffvalues(cfit);
        fit_line = coeffs(1)*cells_ci_loc_y_l*pix2um+coeffs(2);
        plot(cells_ci_loc_y_l*pix2um,fit_line,'m','LineWidth',1.5);                        
        title("Left Boundary")
        xlabel("Distance um")
        ylabel("Basal level")
        xlim([0 400])    
        hold off
    subplot(1,2,2)
        hold on
        scatter(cells_ci_loc_y_r*pix2um,basal_level_r,'c')
        cfit = fit(cells_ci_loc_y_r*pix2um,basal_level_r','poly1');    
        coeffs = coeffvalues(cfit);
        fit_line = coeffs(1)*cells_ci_loc_y_r*pix2um+coeffs(2);
        plot(cells_ci_loc_y_r*pix2um,fit_line,'m','LineWidth',1.5);         
        title("Right Boundary")
        xlabel("Distance um")
        ylabel("Basal level")
        xlim([0 400])    
        hold off
    sgtitle("Basal level Vs Y")


      % figure
      % [X,Y] = meshgrid(cells_ci_d2inj_l,cells_ci_loc_y_l*pix2um);
      % Z=zeros(size(X));
      % surf(X,Y,basal_level_l)

    % surf(cells_ci_d2inj_l,cells_ci_loc_y_l*pix2um,basal_level_l)

    coeffs = [];
    if en_fit_analysis
        n_synth_prof= 5;
        fitmodel = 'secord';
        colorsm = rand(n_synth_prof,3);
        colorsm = [62, 120, 229;36, 230, 45;230, 110, 61;130, 54, 229;79, 98, 102]/255;
        % NCLS(class_i).fitParams = classFitStats(profiles_i(:,frame_start:frame_end)',20,fitmodel);        
        NCLS(class_i).fitmodel = fitmodel;
        coeffs_matrix = genCoeffsMatrix(NCLS(class_i).fitParams,n_synth_prof);
        % [n_coeff,n_params] = size(NCLS(1).fitParams);
        % coeffs_matrix = zeros(n_synth_prof,n_coeff);
        % for i_coeff=1:1:n_coeff
        %     coeff_mu = NCLS(class_i).fitParams(i_coeff,1);
        %     coeff_sd = NCLS(class_i).fitParams(i_coeff,2);
        %     coeffs_matrix(:,i_coeff) = coeff_mu + coeff_sd * randn(n_synth_prof,1); % Generation of coefficients with extracted variability
        % end
        figure       
        hold on
        if ~isempty(lim_global)
            dist = linspace(5,lim_global,n_synth_prof);    
        else
            dist = linspace(10,50,n_synth_prof);
        end        

        for i_sp=1:1:n_synth_prof                       
            % basal ==> basal
            coeffs(1) = genProfileFit(NCLS(class_i).fitparams_basal,dist(i_sp),NCLS(class_i).fitparams_basal_model,0,0);
            % prominence ==> zeta
            coeffs(2) = genProfileFit(NCLS(class_i).fitparams_prominence,dist(i_sp),NCLS(class_i).fitparams_prominence_model,0,0);
            % t_r ==> w_n
            coeffs(3) = genProfileFit(NCLS(class_i).fitparams_peaktime,dist(i_sp),NCLS(class_i).fitparams_peaktime_model,0,0);
            % plateau
            coeffs(4) = genProfileFit(NCLS(class_i).fitparams_plateau,dist(i_sp),NCLS(class_i).fitparams_plateau_model,0,0);
            % Plateau probability
            coeffs(5) = NCLS(class_i).prob_plateau;

            coeffs(6) = NCLS(class_i).th_plateu;

            coeffs(7) = NCLS(class_i).plateau_range_low    ;    
            coeffs(8) = NCLS(class_i).plateau_range_medium ;
            coeffs(9) = NCLS(class_i).plateau_range_high   ;
    
            coeffs(10) = NCLS(class_i).plateau_prob_low    ;    
            coeffs(11) = NCLS(class_i).plateau_prob_medium ;
            coeffs(12) = NCLS(class_i).plateau_prob_high   ;


            if (sum(coeffs<0)>0)
                warning(strcat("Got coefficients <0 for d=",num2str(dist(i_sp)),", value was limited to 0, please chech distance ranges."))
                coeffs
                coeffs(coeffs<0)=0;  % crop negative values
            else
                synt_profile = genProfileFit(coeffs,t,NCLS(class_i).fitmodel,0,0);
                disp(strcat("Dist=",num2str(dist(i_sp))))
                plot(t(1:frame_end-frame_start),synt_profile(1:frame_end-frame_start),'Color',colorsm(i_sp,:),'LineWidth',1.2)            
                % mstring = strcat("N=",num2str(i_sp),"  Dist=",num2str(dist(i_sp)),"  NProf=",num2str(i_sp),"  Basal=",num2str(coeffs(1)),"  PeakProm=",num2str(coeffs(2))," PeakTime=",num2str(coeffs(3)));
                mstring = strcat("D2Inj=",num2str(dist(i_sp)),"um  Basal=",num2str(coeffs(1)),"  MP=",num2str(coeffs(2))," MP_{time}=",num2str(coeffs(3)),"s");
                % text_time = floor(length(t)/3);
                % text(text_time,max(synt_profile),mstring,'Color',colorsm(i_sp,:),'FontSize',8)                
                text_time = 0.12*t(end);
                % text(text_time,synt_profile(end)+0.012,mstring,'Color',colorsm(i_sp,:),'FontSize',14,'Interpreter','tex')
                text(text_time,1.8-(i_sp-1)*0.055,mstring,'Color',colorsm(i_sp,:),'FontSize',14,'Interpreter','tex')
            end 
        end
        hold off
        % title("Class 2: Synthetic Profiles Samples")
        xlabel("Time (seconds)",'Interpreter','tex','FontSize',size_label)
        ylabel("Ratio (F_{340}/F_{380})",'Interpreter','tex','FontSize',size_label)
        xline(50,'--','t_{inj}=50s','FontSize',16,'Interpreter','tex','LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
        set(gca,'FontSize',size_axis)
        ylim([0.78, 1.45]) 
     
        figure       
        hold on
        if ~isempty(lim_global)
            dist = linspace(5,lim_global,n_synth_prof);    
        else
            dist = linspace(10,50,n_synth_prof);
        end        

        for i_sp=1:1:n_synth_prof                       
            % basal ==> basal
            coeffs(1) = genProfileFit(NCLS(class_i).fitparams_basal,dist(i_sp),NCLS(class_i).fitparams_basal_model,0,0);
            % prominence ==> zeta
            coeffs(2) = genProfileFit(NCLS(class_i).fitparams_prominence,dist(i_sp),NCLS(class_i).fitparams_prominence_model,0,0);
            % t_r ==> w_n
            coeffs(3) = genProfileFit(NCLS(class_i).fitparams_peaktime,dist(i_sp),NCLS(class_i).fitparams_peaktime_model,0,0);
            % plateau
            coeffs(4) = genProfileFit(NCLS(class_i).fitparams_plateau,dist(i_sp),NCLS(class_i).fitparams_plateau_model,0,0);
            % Plateau probability
            coeffs(5) = NCLS(class_i).prob_plateau;

            coeffs(6) = NCLS(class_i).th_plateu;

            coeffs(7) = NCLS(class_i).plateau_range_low    ;    
            coeffs(8) = NCLS(class_i).plateau_range_medium ;
            coeffs(9) = NCLS(class_i).plateau_range_high   ;
    
            coeffs(10) = NCLS(class_i).plateau_prob_low    ;    
            coeffs(11) = NCLS(class_i).plateau_prob_medium ;
            coeffs(12) = NCLS(class_i).plateau_prob_high   ;


            if (sum(coeffs<0)>0)
                warning(strcat("Got coefficients <0 for d=",num2str(dist(i_sp)),", value was limited to 0, please chech distance ranges."))
                coeffs
                coeffs(coeffs<0)=0;  % crop negative values
            else
                synt_profile = genProfileFit(coeffs,t,NCLS(class_i).fitmodel,0,0);
                disp(strcat("Dist=",num2str(dist(i_sp))))
                plot(t(1:frame_end-frame_start),synt_profile(1:frame_end-frame_start),'Color',colorsm(i_sp,:),'LineWidth',1.2)            
                % mstring = strcat("N=",num2str(i_sp),"  Dist=",num2str(dist(i_sp)),"  NProf=",num2str(i_sp),"  Basal=",num2str(coeffs(1)),"  PeakProm=",num2str(coeffs(2))," PeakTime=",num2str(coeffs(3)));
                mstring = strcat("D2Inj=",num2str(dist(i_sp)),"um  Basal=",num2str(coeffs(1)),"  MP=",num2str(coeffs(2))," MP_{time}=",num2str(coeffs(3)),"s");
                % text_time = floor(length(t)/3);
                % text(text_time,max(synt_profile),mstring,'Color',colorsm(i_sp,:),'FontSize',8)                
                text_time = 0.12*t(end);
                % text(text_time,synt_profile(end)+0.012,mstring,'Color',colorsm(i_sp,:),'FontSize',14,'Interpreter','tex')
                text(text_time,1.8-(i_sp-1)*0.055,mstring,'Color',colorsm(i_sp,:),'FontSize',14,'Interpreter','tex')
            end 
        end
        hold off
        % title("Class 2: Synthetic Profiles Samples")
        xlabel("Time (seconds)",'Interpreter','tex','FontSize',size_label)
        ylabel("Ratio (F_{340}/F_{380})",'Interpreter','tex','FontSize',size_label)
        xline(50,'--','t_{inj}=50s','FontSize',16,'Interpreter','tex','LabelVerticalAlignment','top','LabelHorizontalAlignment','center')
        set(gca,'FontSize',size_axis)
        ylim([0.78, 1.45])
        % NCLS(class_i).fit_prof_mean  = genProfileFit(NCLS(class_i).fitParams(:,1)',t,fitmodel,0,0);
    end


end
