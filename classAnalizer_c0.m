function [NCLS] = classAnalizer_c0(NCLS,cells_online,ca_prof_smooth,params,inj_band,pix2um,fitmodel)

n_online = max(size(cells_online));
Ts = 3; % TODO: Connect Ts to highest level in case user change this value from GUI
if params.tstart==0
    frame_start = 1;    
else
    frame_start = floor(params.tstart/Ts);
end
frame_inj = floor(params.tstim/Ts);
frame_end = floor(params.tend/Ts);


band_offset =  floor((inj_band(end)-inj_band(1))/2);
en_fit_analysis = 0;

class_number = max(size(NCLS));

f = waitbar(0,'Please wait...');

for class_i = 1    
    cells_class_i = cells_online(NCLS(class_i).cVector);      % Filter only Class 2 cells    
    n_prof_i = max(size(cells_class_i));
    n_cells_i = length(cells_class_i);
    n_frames_i = max(size(cells_class_i(1).ca_profile));            % only valid since Class2 are online cells
    t = (0:1:n_frames_i-1)*Ts;

    profiles_i = reshape([cells_class_i(:).ca_profile],n_frames_i,[])'; 
     
    basal_lines = mean(profiles_i(:,frame_start:frame_inj),2);
    basal_level.mean = mean(basal_lines);
    basal_level.std  = std(basal_lines);
    % 

   %% Computation: basal level per distance
    i_cell_r = 0; i_cell_l = 0;
    cells_ci_loc_x_l = [];cells_ci_loc_x_r = []; 
    cells_ci_d2inj_l = []; cells_ci_d2inj_r = []; cells_ci_d2inj = [];  % X
    cells_ci_loc_y_l = []; cells_ci_loc_y_r = []; % Y
    basal_level_r = []; basal_level_l = [];

    for i_cell = 1:1:n_cells_i
            cells_ci_loc(i_cell,1:2) = cells_class_i(i_cell).xy_hist(1,:);
            if not(isempty(inj_band))
                % cells_ci_loc_y(j,i_cell) = cells_ci_loc(i_cell,2); 
                if inj_band(1)>=cells_ci_loc(i_cell,1)
                    i_cell_l = i_cell_l+1;
                    cells_ci_loc_x_l(i_cell_l,:) = cells_ci_loc(i_cell,:); 
                    cells_ci_loc_y_l(i_cell_l,:) = cells_ci_loc(i_cell,2); 
                    cells_ci_d2inj(i_cell) = abs(inj_band(1)-cells_ci_loc(i_cell,1))*pix2um;
                    basal_level_l(i_cell_l) =  mean(profiles_i(i_cell,frame_start:frame_inj));
                elseif inj_band(end)<=cells_ci_loc(i_cell,1)
                    i_cell_r = i_cell_r+1;
                    cells_ci_loc_x_r(i_cell_r,:) = cells_ci_loc(i_cell,:); 
                    cells_ci_loc_y_r(i_cell_r,:) = cells_ci_loc(i_cell,2); 
                    cells_ci_d2inj(i_cell) = abs(inj_band(end)-cells_ci_loc(i_cell,1))*pix2um;
                    basal_level_r(i_cell_r) =  mean(profiles_i(i_cell,frame_start:frame_inj));                    
                else
                    cells_ci_loc(i_cell,:) = [];
                    cells_ci_d2inj(i_cell) = 9999;
                end
            end
    end
    
    cells_ci_d2inj_l = abs((inj_band(1)+band_offset)-cells_ci_loc_x_l(:,1))*pix2um;    % Left side distance to inj center
    cells_ci_d2inj_r = abs((inj_band(end)-band_offset)-cells_ci_loc_x_r(:,1))*pix2um; % Rigth side distance to inj center

    cells_ci_d2inj_l_b = abs((inj_band(1))-cells_ci_loc_x_l(:,1))*pix2um;    % Left side distance to inj boundary
    cells_ci_d2inj_r_b = abs(cells_ci_loc_x_r(:,1)-(inj_band(end)))*pix2um; % Rigth side distance to inj boundary     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xmax = max([max(cells_ci_d2inj_l),max(cells_ci_d2inj_r)])+10;
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
    figure
    subplot(1,3,1)
        hold on
        scatter(cells_ci_d2inj_l,basal_level_l,'c')        
        rsquare = [];                   % 2 report
        rms = [];                       % 2 report
        class_feature = 'C0_Basal_L';   % 2 report
        if en_several_models == 1                         
            for i_model = 1:1:n_models
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_l_b,basal_level_l',fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report
                rms(i_model,1) = gof.rmse;          % 2 report
                cl_ft{i_model,:} = class_feature;   % 2 report
                coeffs1 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs1,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            legend(model_list)

            cfit = fit(cells_ci_d2inj_l_b,basal_level_l',fitmodel_basal);
            coeffs_left = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs_left,cells_ci_d2inj_l_b,fitmodel_basal,0,0);
            plot(cells_ci_d2inj_l_b,fit_line_l,'m','LineWidth',1.5);
        else
            cfit = fit(cells_ci_d2inj_l_b,basal_level_l',fitmodel_basal);
            coeffs_left = coeffvalues(cfit);
            fit_line_l = genProfileFit(coeffs_left,cells_ci_d2inj_l_b,fitmodel_basal,0,0);
            plot(cells_ci_d2inj_l_b,fit_line_l,'m','LineWidth',1.5);
        end
        % cfit = fit(cells_ci_d2inj_l,basal_level_l',fitmodel_basal);    
        % coeffs_left = coeffvalues(cfit);
        % % fit_line_l = coeffs_left(1)*cells_ci_d2inj_l+coeffs_left(2);
        % fit_line_l = genProfileFit(coeffs_left,cells_ci_d2inj_l,fitmodel_basal,0,0);
        % plot(cells_ci_d2inj_l,fit_line_l,'m',LineWidth=1.5);        
        title("Left Boundary")
        xlabel("Distance um")
        ylabel("Basal level")
        xlim([0 400])
        ylim([0.5 1.5])
        hold off
    subplot(1,3,2)
        hold on
        scatter(cells_ci_d2inj,basal_lines,'c')
        cfit = fit(cells_ci_d2inj',basal_lines,fitmodel_basal);    
        coeffs = coeffvalues(cfit);       
        % fit_line = coeffs(1)*cells_ci_d2inj+coeffs(2);
        fit_line = genProfileFit(coeffs,cells_ci_d2inj,fitmodel_basal,0,0);        
        NCLS(class_i).fitbasal_params = coeffs;
        plot(cells_ci_d2inj,fit_line,'m',LineWidth=1.5);        
        title("Center")
        xlabel("Distance um")
        ylabel("Basal level")
        xlim([0 400])
        ylim([0.5 1.5])
        hold off
    subplot(1,3,3)
        hold on
        scatter(cells_ci_d2inj_r_b,basal_level_r,'c')
        rsquare = [];                   % 2 report
        rms = [];                       % 2 report
        class_feature = 'C0_Basal_R';% 2 report        
        model_list(1)=[];
          if en_several_models == 1 
            for i_model = 1:1:n_models
                fmodel = model_list{i_model};
                [cfit, gof] = fit(cells_ci_d2inj_r_b,basal_level_r',fmodel);
                rsquare(i_model,1) = gof.rsquare;   % 2 report
                rms(i_model,1) = gof.rmse;          % 2 report
                cl_ft{i_model,:} = class_feature;   % 2 report
                coeffs2 = coeffvalues(cfit);
                fit_line = genProfileFit(coeffs2,distv,fmodel,0,0);
                plot(distv,fit_line,'LineWidth',1.5,'Color',colorsm(i_model,:))
            end
            tableData = table(cl_ft,model_list',rsquare,rms,'VariableNames',{'ClassFeature','Model','R2','RMS'});    % 2 report
            writetable(tableData,'fitReport.txt','Delimiter','\t','WriteMode','append');                             % 2 report            
            model_list{end+1} = 'ExpData';
            model_list = circshift(model_list,1);
            model_list{end+1} = 'Mean Fit';
            legend(model_list)
            
            cfit = fit(cells_ci_d2inj_r_b,basal_level_r',fitmodel_basal);
            coeffs_rigth = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs_rigth,cells_ci_d2inj_r_b,fitmodel_basal,0,0);
            plot(cells_ci_d2inj_r_b,fit_line_r,'m','LineWidth',1.5);
        else
            cfit = fit(cells_ci_d2inj_r_b,basal_level_r',fitmodel_basal);
            coeffs_rigth = coeffvalues(cfit);
            fit_line_r = genProfileFit(coeffs_rigth,cells_ci_d2inj_r_b,fitmodel_basal,0,0);
            plot(cells_ci_d2inj_r_b,fit_line_r,'m','LineWidth',1.5);
          end
        % cfit = fit(cells_ci_d2inj_r,basal_level_r',fitmodel_basal);    
        % coeffs_rigth = coeffvalues(cfit);
        % % fit_line_r = coeffs_rigth(1)*cells_ci_d2inj_r+coeffs_rigth(2);
        % fit_line_r = genProfileFit(coeffs_rigth,cells_ci_d2inj_r,fitmodel_basal,0,0);
        % plot(cells_ci_d2inj_r,fit_line_r,'m',LineWidth=1.5)                
        title("Right Boundary")
        xlabel("Distance um")
        ylabel("Basal level")
        xlim([0 400])
        ylim([0.5 1.5])
        hold off
    sgtitle("Basal level Vs Distance to injury")



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
    plot(cells_ci_d2inj_l_b,fit_line_l,'--',LineWidth=1.5,Color=color_l);        
    scatter(cells_ci_d2inj_r_b,basal_level_r,55,'MarkerEdgeColor',color_r,'Marker',marker_r) % Rigth data
    plot(cells_ci_d2inj_r_b,fit_line_r,'--',LineWidth=1.5,Color=color_r);
    coeffs_mean = mean([coeffs_left;coeffs_rigth],1);   % 2 model
    fit_line_m_basal = genProfileFit(coeffs_mean,distv,fitmodel_basal,0,0);
    plot(distv,fit_line_m_basal,'-',LineWidth=1.5,Color=color_m)
      NCLS(class_i).fitparams_basal_model = fitmodel_basal;        
      NCLS(class_i).fitbasal_params = coeffs_mean;
    xlabel("Distance to injury (\mum)",'Interpreter','tex','FontSize',size_label)
    ylabel("Basal level (F_{340}/F_{380})",'Interpreter','tex','FontSize',size_label)
    set(gca,'FontSize',size_axis)
    xlim([0 max([max(cells_ci_d2inj_l_b),max(cells_ci_d2inj_r_b)])+10])
    % ylim([min([basal_level_l,basal_level_r])-0.1 max([basal_level_l,basal_level_r])+0.1])
    ylim([0.7 1.15]);
    qw{1} = plot(nan,"o",'MarkerEdgeColor',color_l);
    qw{2} = plot(nan,"+",'MarkerEdgeColor',color_r);
    qw{3} = plot(nan, '-',LineWidth=1.5,Color=color_m);
    legend([qw{:}],{'Left side data','Rigth side data','Model'},'FontSize',size_legend)    
    hold off

    FeatureName = "C0_BasalLevel";                                                                                 % 2 report B
    tableData = table(FeatureName,{fitmodel_basal},coeffs_mean,'VariableNames',{'Feature','Model','Coeffs'});        % 2 report B
    writetable(tableData,'fitModels.txt','Delimiter','\t','WriteMode','append');                                   % 2 report B
    



    figure
    subplot(1,2,1)
        hold on
        scatter(cells_ci_loc_y_l*pix2um,basal_level_l,'c')
        cfit = fit(cells_ci_loc_y_l*pix2um,basal_level_l','poly1');    
        coeffs = coeffvalues(cfit);
        fit_line = coeffs(1)*cells_ci_loc_y_l*pix2um+coeffs(2);
        plot(cells_ci_loc_y_l*pix2um,fit_line,'m',LineWidth=1.5);                        
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
        plot(cells_ci_loc_y_r*pix2um,fit_line,'m',LineWidth=1.5);         
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

        NCLS(class_i).fitParams = classFitStats(profiles_i',20,fitmodel);
        NCLS(class_i).fitmodel = fitmodel;

    if en_fit_analysis
        n_synth_prof= 30;
        NCLS(class_i).fitParams = classFitStats(profiles_i',20,fitmodel);
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
        for i_sp=1:1:n_synth_prof            
            synt_profile = genProfileFit(coeffs_matrix(i_sp,:),t,fitmodel,0,0);
            plot(t(1:frame_end-frame_start),synt_profile(1:frame_end-frame_start))                        
        end
        hold off
        title("Class 0: Synthetic Profiles Samples")        
        NCLS(class_i).fit_prof_mean  = genProfileFit(NCLS(class_i).fitParams(:,1)',t,fitmodel,0,0);
    end

    bar_message = strcat("Extracting stats Class 0 ",num2str(class_i));
    waitbar(class_i/class_number,f,bar_message);
end

close(f)