function sytheticFramePlotSeq(sim_cell,sim_canv,video_en,video_path,inj_x)
thesis_plot_en = 1;
if video_en
    
    theFiles = dir(fullfile(video_path,'images'));
    for k = 1 : length(theFiles)
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(fullfile(video_path,'images'), baseFileName);
        fprintf(1, 'Now deleting %s\n', fullFileName);
        delete(fullFileName);
    end
     
    workingDir = video_path;
    mkdir(workingDir)
    mkdir(workingDir,'images')
end
[m_canv,~]=size(sim_canv)
n_scells = max(size(sim_cell));
n_frames = length([sim_cell(1).icp])

%% Thesis plot begin
if thesis_plot_en
    size_axis = 20;
    size_label = 24;
    size_legend = 14;
    color_m = [15,15,15]/255;
    
    colors = [0.08,0.08,0.08;
        0.3,0.3,0.3;
        0.6,0.6,0.6;
        0.85,0.85,0.85];
        
    figure;
    hold on
    for k=1:1:n_scells   
             if ((sim_cell(k).class == 4 ) ||(sim_cell(k).class == 3 ) ||(sim_cell(k).class == 2 ) || (sim_cell(k).class == 1 ))
             plot([1:1:n_frames],sim_cell(k).icp,'LineWidth',1,'Color',colors(sim_cell(k).class,:))
             xlabel("Frame",'FontSize',size_label)
             ylabel("Normalized Intensity",'FontSize',size_label)
             %set(gca,'FontSize',size_axis)        
             end
    end
    xlim([1 n_frames])
    set(gca,'FontSize',size_axis)
    qw{1} = plot(nan,"-",LineWidth=1,Color=colors(1,:));
    qw{2} = plot(nan,"-",LineWidth=1,Color=colors(2,:));
    qw{3} = plot(nan,"-",LineWidth=1,Color=colors(3,:));
    qw{4} = plot(nan,"-",LineWidth=1,Color=colors(4,:));
    legend([qw{:}],{'Class0','Class1','Class2','Class3'},'FontSize',size_legend)    
    hold off
end
%% Thesis plot end

c0v = [];
c1v = [];
c2v = [];
c3v = [];
% Rescale intensity to range [0,1]
icp_max = max([sim_cell(:).icp]);
icp_min = min([sim_cell(:).icp]);    
for k=1:1:n_scells
    sim_cell(k).icp = (sim_cell(k).icp-icp_min)/(icp_max-icp_min);

    if sim_cell(k).class == 1
        disp('gotc0')
        c0v=[c0v k];
    end
    if sim_cell(k).class == 2
        disp('gotc1')
        c1v=[c1v k];
    end
    if sim_cell(k).class == 3
        disp('gotc2')
        c2v=[c2v k];
    end    
    if sim_cell(k).class == 4
        disp('gotc3')
        c3v=[c3v k];
    end
    % if sim_cell(k).class == 3
    %     figure;
    %     plot([1:1:n_frames],sim_cell(k).icp)
    %     xlabel("Frame")
    %     ylabel("Ratio")
    %     title("Class2: Synthetic profile - syntheticFramePlotSeq")
    % end
end

%% Thesis plot begin
if thesis_plot_en
foi = [23,24,25,26,27,28,29,30,31,32];

if ~isempty(c0v)
    c0_idx = c0v(round(length(c0v)/2));
    disp('c0_location')
    sim_cell(c0_idx).params
    for i_foi = 1:1:length(foi)
        figure;  
        hold on
        plot([1:1:n_frames],sim_cell(c0_idx).icp,'LineWidth',2,'Color',colors(1,:));
        plot(foi,sim_cell(c0_idx).icp(foi),'k','Marker','v','MarkerSize',8);
        plot(foi(i_foi),sim_cell(c0_idx).icp(foi(i_foi)),'r','Marker','v','MarkerSize',12,'MarkerFaceColor','r');
        hold off
        xlim([foi(1)-5 foi(end)+5])
        ylim([0 1])
        xlabel("Frame",'FontSize',size_label)
        ylabel("Normalized Intensity",'FontSize',size_label)
        set(gca,'FontSize',size_axis)
    end
    figure;
    hold on
    plot([1:1:n_frames],sim_cell(c0_idx).icp,'LineWidth',1,'Color',colors(1,:));
    plot(foi,sim_cell(c0_idx).icp(foi),'LineWidth',2.5,'Color',[1 0 0]);    
    hold off
    xlim([1 n_frames])
    ylim([0 1])
    xlabel("Frame",'FontSize',size_label)
    ylabel("Normalized Intensity",'FontSize',size_label)
    set(gca,'FontSize',size_axis)       
end

if ~isempty(c1v)
    c1_idx = c1v(round(length(c1v)/2));
    disp('c1_location')
    sim_cell(c1_idx).params
    for i_foi = 1:1:length(foi)
        figure;  
        hold on
        plot([1:1:n_frames],sim_cell(c1_idx).icp,'LineWidth',2,'Color',colors(1,:));
        plot(foi,sim_cell(c1_idx).icp(foi),'k','Marker','v','MarkerSize',8);
        plot(foi(i_foi),sim_cell(c1_idx).icp(foi(i_foi)),'r','Marker','v','MarkerSize',12,'MarkerFaceColor','r');        
        hold off
        xlim([foi(1)-5 foi(end)+5])
        ylim([0 1])
        xlabel("Frame",'FontSize',size_label)
        ylabel("Normalized Intensity",'FontSize',size_label)
        set(gca,'FontSize',size_axis)        
    end 
    figure;
    hold on
    plot([1:1:n_frames],sim_cell(c1_idx).icp,'LineWidth',1,'Color',colors(1,:));
    plot(foi,sim_cell(c1_idx).icp(foi),'LineWidth',2,'Color',[1 0 0]);    
    hold off
    xlim([1 n_frames])
    ylim([0 1])
    xlabel("Frame",'FontSize',size_label)
    ylabel("Normalized Intensity",'FontSize',size_label)
    set(gca,'FontSize',size_axis)   
end

if ~isempty(c2v)
    c2_idx = c2v(round(length(c2v)/2));
    disp('c2_location')
    sim_cell(c2_idx).params
    for i_foi = 1:1:length(foi)
        figure;  
        hold on
        plot([1:1:n_frames],sim_cell(c2_idx).icp,'LineWidth',2,'Color',colors(1,:));
        plot(foi,sim_cell(c2_idx).icp(foi),'k','Marker','v','MarkerSize',8);
        plot(foi(i_foi),sim_cell(c2_idx).icp(foi(i_foi)),'r','Marker','v','MarkerSize',12,'MarkerFaceColor','r');
        hold off
        xlim([foi(1)-5 foi(end)+5])
        ylim([0 1])
        xlabel("Frame",'FontSize',size_label)
        ylabel("Normalized Intensity",'FontSize',size_label)
        set(gca,'FontSize',size_axis)        
    end
    figure;
    hold on
    plot([1:1:n_frames],sim_cell(c2_idx).icp,'LineWidth',1,'Color',colors(1,:));
    plot(foi,sim_cell(c2_idx).icp(foi),'LineWidth',2,'Color',[1 0 0]);    
    hold off
    xlim([1 n_frames])
    ylim([0 1])
    xlabel("Frame",'FontSize',size_label)
    ylabel("Normalized Intensity",'FontSize',size_label)
    set(gca,'FontSize',size_axis)    
end

if ~isempty(c3v)
    c3_idx = c3v(round(length(c3v)/2));
    disp('c3_location')
    sim_cell(c3_idx).params
    for i_foi = 1:1:length(foi)
        figure;  
        hold on
        plot([1:1:n_frames],sim_cell(c3_idx).icp,'LineWidth',2,'Color',colors(1,:));
        plot(foi,sim_cell(c3_idx).icp(foi),'k','Marker','v','MarkerSize',8);
        plot(foi(i_foi),sim_cell(c3_idx).icp(foi(i_foi)),'r','Marker','v','MarkerSize',12,'MarkerFaceColor','r');
        hold off
        xlim([foi(1)-5 foi(end)+5])
        ylim([0 1])
        xlabel("Frame",'FontSize',size_label)
        ylabel("Normalized Intensity",'FontSize',size_label)
        set(gca,'FontSize',size_axis)        
    end
    figure;
    hold on
    plot([1:1:n_frames],sim_cell(c3_idx).icp,'LineWidth',1,'Color',colors(1,:));
    plot(foi,sim_cell(c3_idx).icp(foi),'LineWidth',2,'Color',[1 0 0]);    
    hold off
    xlim([1 n_frames])
    ylim([0 1])
    xlabel("Frame",'FontSize',size_label)
    ylabel("Normalized Intensity",'FontSize',size_label)
    set(gca,'FontSize',size_axis)    
end

end
%% Thesis plot end

%% Thesis plot begin
if thesis_plot_en
    size_axis = 20;
    size_label = 24;
    size_legend = 14;
    color_m = [15,15,15]/255;
    
    colors = [0.08,0.08,0.08;
        0.3,0.3,0.3;
        0.6,0.6,0.6;
        0.85,0.85,0.85];
        
    figure;
    hold on
    for k=1:1:n_scells   
             if ((sim_cell(k).class == 4 ) ||(sim_cell(k).class == 3 ) ||(sim_cell(k).class == 2 ) || (sim_cell(k).class == 1 ))
             plot([1:1:n_frames],sim_cell(k).icp,'LineWidth',1,'Color',colors(sim_cell(k).class,:))
             xlabel("Frame",'FontSize',size_label)
             ylabel("Ratio (F_{340}/F_{380})",'FontSize',size_label)
             %set(gca,'FontSize',size_axis)        
             end
    end
    xlim([1 n_frames])
    set(gca,'FontSize',size_axis)
    qw{1} = plot(nan,"-",LineWidth=1,Color=colors(1,:));
    qw{2} = plot(nan,"-",LineWidth=1,Color=colors(2,:));
    qw{3} = plot(nan,"-",LineWidth=1,Color=colors(3,:));
    qw{4} = plot(nan,"-",LineWidth=1,Color=colors(4,:));
    legend([qw{:}],{'Class0','Class1','Class2','Class3'},'FontSize',size_legend)    
    hold off
end
%% Thesis plot end

%% Thesis plot begin
if thesis_plot_en
    size_axis = 20;
    size_label = 24;
    size_legend = 14;
    color_m = [15,15,15]/255;
    
    colors = [0.08,0.08,0.08;
        0.3,0.3,0.3;
        0.6,0.6,0.6;
        0.85,0.85,0.85];
        
    figure;
    hold on
    for k=1:1:25   
             if ((sim_cell(k).class == 4 ) ||(sim_cell(k).class == 3 ) ||(sim_cell(k).class == 2 ) || (sim_cell(k).class == 1 ))
             plot([1:1:n_frames],sim_cell(k).icp,'LineWidth',1,'Color',colors(sim_cell(k).class,:))
             xlabel("Frame",'FontSize',size_label)
             ylabel("Ratio (F_{340}/F_{380})",'FontSize',size_label)
             %set(gca,'FontSize',size_axis)        
             end
    end
    xlim([1 n_frames])
    set(gca,'FontSize',size_axis)
    qw{1} = plot(nan,"-",LineWidth=1,Color=colors(1,:));
    qw{2} = plot(nan,"-",LineWidth=1,Color=colors(2,:));
    qw{3} = plot(nan,"-",LineWidth=1,Color=colors(3,:));
    qw{4} = plot(nan,"-",LineWidth=1,Color=colors(4,:));
    legend([qw{:}],{'Class0','Class1','Class2','Class3'},'FontSize',size_legend)    
    hold off
end
%% Thesis plot end


figure
for kf=1:n_frames
    % figure;
    imshow(sim_canv);
    if kf == 1
        set(gca,'position',[0 0 1 1],'units','normalized')
%         ax = gca;
%         outerpos = ax.OuterPosition;
%         ti = ax.TightInset; 
%         left = outerpos(1) + ti(1);
%         bottom = outerpos(2) + ti(2);
%         ax_width = outerpos(3) - ti(1) - ti(3);
%         ax_height = outerpos(4) - ti(2) - ti(4);
%         ax.Position = [left bottom ax_width ax_height];
    end    

    hold on;
    for k=1:1:n_scells        
        % rectangle('Position',sim_cell(k).params,'Curvature',[1 1],'EdgeColor',[sim_cell(k).icp(kf) 0 0],"FaceColor",[sim_cell(k).icp(kf) 0 0]) % Draw s-cell> per cell only for PAPER          
        color_rgb = sim_cell(k).color;
        color_hsv = rgb2hsv(color_rgb);
        color_hsv(3) = sim_cell(k).icp(kf);
        color_rgb = hsv2rgb(color_hsv);
        rectangle('Position',sim_cell(k).params,'Curvature',[1 1],'EdgeColor',sim_cell(k).color,"FaceColor",sim_cell(k).color) % Draw s-cell> per cell   
        rectangle('Position',[sim_cell(k).params(1)+0.125*sim_cell(k).params(3) sim_cell(k).params(2)+0.125*sim_cell(k).params(4) 0.75*sim_cell(k).params(3) 0.75*sim_cell(k).params(4)],'Curvature',[1 1],'EdgeColor',[0 0 0],'LineWidth',0.5,"FaceColor",color_rgb) % Draw s-cell> per cell   
        rectangle('Position',[inj_x 1 1 m_canv-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3) %% injury
    end
    flabel=strcat("f",num2str(kf));    
    text(30,30,flabel,"Color",'g','FontSize',30,'FontWeight','bold');
    hold off;
    if video_en
        filename = [sprintf('%03d',kf)];
        fullname = fullfile(workingDir,'images',filename);
        saveas(gcf,fullname,'jpg')
%         imwrite(img,fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
    end
%     pause(0.01)
end
if video_en
    display("Start Video Gen ...")
    imageNames = dir(fullfile(workingDir,'images','*.jpg'));
    imageNames = {imageNames.name}';
    outputVideo = VideoWriter(fullfile(workingDir,'synthcells5.avi'));
    outputVideo.FrameRate = 20;
    open(outputVideo)
    for ii = 1:length(imageNames)
        img = imread(fullfile(workingDir,'images',imageNames{ii}));
        writeVideo(outputVideo,img)
    end
    close(outputVideo)    
    disp("End Video Gen ...stored in ")
    fullfile(workingDir,'synthcells5.avi')
end