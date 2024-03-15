clear all; close all; clc

path_data_feat = "C:\Users\marci\OneDrive\Documentos\MATLAB\sandbox\Ca_App_V2-main\common_twin_data";
path_data_full = "C:\Users\marci\OneDrive\Documentos\MATLAB\sandbox\Ca_App_V2-main\common_twin_data";
path_results =  "C:\Users\marci\OneDrive\Documentos\MATLAB\sandbox\Ca_App_V2-main\common_twin_data\StatAnalysis";

ExpFull = [      "100528_1";      "100528_4";      "100528_5";      "100528_6";      "100528_8";      "100601_1";];
ExpFeat = ["CFeat_100528_1";"CFeat_100528_4";"CFeat_100528_5";"CFeat_100528_6";"CFeat_100528_8";"CFeat_100601_1";];

n_exp = length(ExpFeat);
subplot = @(m,n,p) subtightplot (m, n, p, [0.055 0.035], [0.08 0.11], [0.025 0.02]);

glb_cells_ci_d_y_um = []; hs_glb_ci_x = [];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   LOOP START  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_exp = 1:1:n_exp
    % clearvars -except i_exp path_data path_results Exp gf_cells_x gf_cells_y gf_cluster_class_data
    % close all; clc
    %% %%%%%%%%%%%%%%%%   CLEAN UP  %%%%%%%%%%%%%%%%%%%
    disp("Clean up start")    
        hstack_ci_y = []; hstack_norm_ci_y=[]; hist_stack = [];
    disp("Clean up end")    
    %% %%%%%%%%%%%%%%%%   LOAD DATA  %%%%%%%%%%%%%%%%%%%    
    disp(strcat("Dislpaying: ",ExpFeat(i_exp)))
    load(fullfile(path_data_feat,ExpFeat(i_exp))) % Loading experiment stats
    load(fullfile(path_data_full,ExpFull(i_exp))) % Loading experiment full data
    %% %%%%%%%%%%%%%%%%   PLOTS: Per Exp %%%%%%%%%%%%%%%
    
    %============ Absolute Ind: X ============%
    if i_exp==1; fig_x = figure; end   
    figure(fig_x)
    subplot(2,n_exp/2,i_exp)
    bar(v_edges(1:end-1),hist_stack,'stacked');xlabel("Distance [um]");ylabel("Freq")
    title(ExpFeat(i_exp))
    sgtitle("X Distributions Normalized %")
    %--------- Stacked Global data collection : X
    if i_exp == 1
        hs_glb_ci_x = hist_stack; 
    else
        hs_glb_ci_x = hs_glb_ci_x+hist_stack;        
    end    

    %============ Absolute Ind: Y ============%
    if i_exp==1; fig_y = figure; end
    yn_classes = length(NCLS);
    nbins = 20;
    v_edges_y = linspace(0,600,2*nbins);
    v_edges_um_y = v_edges_y*pix2um;
    cells_ci_d_y_um = cells_ci_d_y.*pix2um;        
    for i_class = 1:1:yn_classes
        figure(fig_y)
        subplot(n_exp,yn_classes,yn_classes*(i_exp-1)+i_class);
        % histogram(nonzeros(cells_ci_d_y(i_class,:)),'BinEdges',v_edges,"FaceColor",NCLS(i_class).color/255); ylabel("Frequency"); xlabel("Distance [um]"); %pixels 
        histogram(nonzeros(cells_ci_d_y_um(i_class,:)),'BinEdges',v_edges_um_y,"FaceColor",NCLS(i_class).color/255); ylabel("Frequency"); xlabel("Distance [um]"); %um
        if i_exp==1
            string_title = strcat("Class ",num2str(i_class));
        else
            string_title = "";
        end
        if i_class == 1 
             string_title = strcat(ExpFeat(i_exp)," - ",string_title);
        end
        title(string_title,'Interpreter','none');        
        xlim([0,600*pix2um]);
        sgtitle("Y distributions")
        % Stack Ind data collection : Y
        hstack_ci_y(:,i_class)=histcounts(nonzeros(cells_ci_d_y_um(i_class,:)),'BinEdges',v_edges_um_y);
        hstack_norm_ci_y = hstack_ci_y./sum(hstack_ci_y,2)*100;
        hstack_norm_ci_y(isnan(hist_stack_normalized))=0;        
    end
    %  Stack Normalized Ind: Y
    % if i_exp==1; fig_y_stack = figure; end
    % figure(fig_y_stack);
    figure
    % subplot(2,n_exp/2,i_exp)
    barnorm = barh(v_edges_um_y(1:end-1),hstack_norm_ci_y,'stacked');xlabel("Distance [um]");ylabel("Freq")
    xlim([0 110]);
    set(gca, 'YDir','reverse')
    legend('Class_0','Class_1','Class_2','Class_3','Interpreter','none')
    string_title = strcat(ExpFull(i_exp)," - ",string_title);
    title(strcat("Y: ",string_title),'Interpreter','none');
    for i_bar =1:1:yn_classes
        barnorm(i_bar).FaceColor = NCLS(i_bar).color./255;
    end
    % sgtitle("Y Distributions Normalized %")
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                  Global Data Collection
    glb_cells_ci_d_y_um = [glb_cells_ci_d_y_um,cells_ci_d_y_um];    % global Y 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                  Cluster Features
    for i_class = 1:1:yn_classes
        plot_features(class_data(i_class),NCLS(i_class).color/255,NCLS(i_class).title)
    end
    disp("pause...")
    pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%   LOOP END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%            GLOBAL PLOTS                   %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% ============    Absolute Global: Y
    fig_glb_y = figure;  
    hst_glb_y = [];
    for i_class = 1:1:yn_classes
        figure(fig_glb_y)        
        subplot(1,yn_classes,i_class);        
        histogram(nonzeros(glb_cells_ci_d_y_um(i_class,:)),'BinEdges',v_edges_um_y,"FaceColor",NCLS(i_class).color/255); ylabel("Frequency"); xlabel("Distance [um]"); 
        string_title = strcat("Class ",num2str(i_class))

        title(string_title,'Interpreter','none');        
        xlim([0,600*pix2um]);
        ylim([0,30]);
        sgtitle("Global: Y distributions")
        
        % Stack Global data collection : Y
        hst_glb_y(:,i_class) = histcounts(nonzeros(glb_cells_ci_d_y_um(i_class,:)),'BinEdges',v_edges_um_y);
        hst_glb_y_norm = hst_glb_y./sum(hst_glb_y,2)*100;
        hst_glb_y_norm(isnan(hst_glb_y_norm))=0;        
    end 
     % Stack Normalized Global: Y
     fig_glb_y_stack = figure;
     barnorm = barh(v_edges_um_y(1:end-1),hst_glb_y_norm,'stacked');xlabel("Distance [um]");ylabel("Freq")
     legend('Class_0','Class_1','Class_2','Class_3','Interpreter','none')
     for i_bar =1:1:yn_classes
         barnorm(i_bar).FaceColor = NCLS(i_bar).color./255;
     end
     sgtitle("Y Global Distributions Normalized %")

     %% ============   Stack Global: X
     for i_class = 1:1:yn_classes      
        % Stack Global data collection : X
        % hst_glb_x(:,i_class) = histcounts(nonzeros(hs_glb_ci_x(:,i_class)),'BinEdges',v_edges);        
        hst_glb_x_norm(:,i_class) = hs_glb_ci_x(:,i_class)./sum(hs_glb_ci_x,2)*100;
        hst_glb_x_norm(isnan(hst_glb_x_norm))=0;                
     end 
     fig_glb_x_stack = figure;
     barnorm = bar(v_edges(1:end-1),hst_glb_x_norm,'stacked');xlabel("Distance [um]");ylabel("Freq")
     legend('Class_0','Class_1','Class_2','Class_3','Interpreter','none')
     for i_bar =1:1:yn_classes
         barnorm(i_bar).FaceColor = NCLS(i_bar).color./255;
     end
     sgtitle("X Global Distributions Normalized %")     

     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

