close all;clear all;clc
visual_paper = 1

video_en=0;video_pathv="C:\Users\marci\OneDrive\Documentos\MATLAB\synthimage";
grid_en = 0;
prompt = {'Rows of cells:','Columns of cells:'};
dlgtitle = 'Matrix size';
dims = [1 35];
definput = {'15','20'};
if visual_paper
    sim_row = 1 
    sim_col = 4 
else
    msize = inputdlg(prompt,dlgtitle,dims,definput);
    sim_row = floor(str2double(msize(1)));
    sim_col = floor(str2double(msize(2)));
end
[sim_cell,sim_canv] = sytheticFrameGen(sim_row,sim_col,grid_en);     % Generate simulation cells & canvas
synt_fig = figure;
figure(synt_fig)
% sytheticFramePlot(sim_cell,sim_canv)
% [inj_x,~] = ginput(1);    % Get injury and plot
% inj_x = round(inj_x);
% [m_canv,~]=size(sim_canv);
% color_matrix = reshape([NCLS(:).color],3,[])'/255;
%     [sim_cell] = syntheticClassAssign(sim_cell,inj_x,color_matrix);
% [sim_cell] = simClassify(sim_cell,inj_x,fis,color_matrix);
[sim_cell] = applyCellModel(sim_cell,"random",sim_row,sim_col);
sytheticFramePlotSeq(sim_cell,sim_canv,video_en,video_pathv)
% rectangle('Position',[inj_x 1 1 m_canv-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3) 

if visual_paper 
    fsize_xylabel = 28;
    fsize_axis = 28;
    lw = 2.5;
    xlim_last = 300 %629
    f1=figure; axes1 = axes('Parent',f1);hold(axes1,'on'); 
    plot([sim_cell(1).icp],'LineWidth',lw);xlim([0 xlim_last])
    set(axes1,'FontName','Arial','FontSize',fsize_axis);
    ylabel('Normalized Intensity','FontSize',fsize_xylabel,'FontName','Arial');
    xlabel('Frame','FontSize',fsize_xylabel,'FontName','Arial');
    hold(axes1,'off');    

    f2=figure; axes2 = axes('Parent',f2);hold(axes2,'on');
    plot([sim_cell(2).icp],'LineWidth',lw);xlim([0 xlim_last])
    set(axes2,'FontName','Arial','FontSize',fsize_axis);
    ylabel('Normalized Intensity','FontSize',fsize_xylabel,'FontName','Arial');
    xlabel('Frame','FontSize',fsize_xylabel,'FontName','Arial');
    hold(axes2,'off');    

    f3=figure; axes3 = axes('Parent',f3);hold(axes3,'on');
    plot([sim_cell(3).icp],'LineWidth',lw);xlim([0 xlim_last])
    set(axes3,'FontName','Arial','FontSize',fsize_axis);
    ylabel('Normalized Intensity','FontSize',fsize_xylabel,'FontName','Arial');
    xlabel('Frame','FontSize',fsize_xylabel,'FontName','Arial');
    hold(axes3,'off');    

    f4=figure; axes4 = axes('Parent',f4);hold(axes4,'on');
    plot([sim_cell(4).icp],'LineWidth',lw);xlim([xlim_last-10 xlim_last])
    set(axes4,'FontName','Arial','FontSize',fsize_axis);
    ylabel('Normalized Intensity','FontSize',fsize_xylabel,'FontName','Arial');
    xlabel('Frame','FontSize',fsize_xylabel,'FontName','Arial');
    hold(axes4,'off');    
end