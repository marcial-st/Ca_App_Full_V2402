clear all;clc;close all;
%% MAIN IDX
i=1;
%% clean report path
% spath = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\ReportGen";
spath = "C:\Users\MST\Documents\ReportGen"
cleanSpath(spath)
%% Main 
cindex = [0,1,2,3,4];
clabels = ["unknown"; %0
    "no response"; %1
    "basal"; %2
    "response"; %3
    "basal+response";%4
    ];

clabels = lower(clabels);
cdict = containers.Map(cindex,clabels);
wfilter= 11;
params =[18,40,0.01,0.25,16,15,wfilter,5,2.75,1,12; %1
         30,40,0.01,0.25,15,15,wfilter,7,2.25,15,25; %2 1.4 std ratio threshold
         50,40,0.01,0.25,10,15,wfilter,7,2.75,15,40; %3 
         40,40,0.01,0.25,10,15,wfilter,7,2.75,15,35; %4 
         65,40,0.01,0.25,25,15,wfilter,7,2.75,20,40; %5 bueno
         65,40,0.01,0.25,25,15,wfilter,7,2.75,10,55; %6 
         60,40,0.01,0.25,15,15,wfilter,7,2.75,5,55; %7 bueno
         45,40,0.01,0.25,10,15,wfilter,7,2.75,10,40; %8 muy pocas muestras, análisis no viable
         60,40,0.01,0.25,10,15,wfilter,7,2.75,20,50]; %9 
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
main_path = "C:\Users\MST\OneDrive\Perfiles\PeakGT";

% gta = ["100526_5_Ajelet.xlsx",...
%     "100528_1_Ajelet.xlsx",...
%     "100528_4_Ajelet.xlsx",...
%     "100528_5_Ajelet.xlsx",...
%     "100528_6_Ajelet.xlsx",...
%     "100528_8_Ajelet.xlsx",...
%     "100601_1_Ajelet.xlsx",...
%     "100604_1_Ajelet.xlsx",...
%     "100604_9_Ajelet.xlsx"];

gta = ["100526_5_Roberto.xlsx",...
    "100528_1_Roberto.xlsx",...
    "100528_4_Roberto.xlsx",...
    "100528_5_Roberto.xlsx",...
    "100528_6_Roberto.xlsx",...
    "100528_8_Roberto.xlsx",...
    "100601_1_Roberto.xlsx",...
    "100604_1_Roberto.xlsx",...
    "100604_9_Roberto.xlsx"];
   
ata = ["ds_100526_5_r6.mat",...
    "ds_100528_1_r6.mat",...
    "ds_100528_4_r6.mat",...
    "ds_100528_5_r6.mat",...
    "ds_100528_6_r6.mat",...
    "ds_100528_8_r6.mat",...
    "ds_100601_1_r6.mat",...
    "ds_100604_4_r6.mat",...    
    "ds_100604_9_r6.mat"];

% gt_path = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Peak_Analysis_Ground_Truth";
% gt_peaks = ["data_100526_5_Roberto.mat",...
%     "data_100528_1_Roberto.mat",...
%     "data_100528_4_Roberto.mat",...
%     "data_100528_5_Roberto.mat",...
%     "data_100528_6_Roberto.mat",...
%     "data_100528_8_Roberto.mat",...
%     "data_100601_1_Roberto.mat",...
%     "data_100604_1_Roberto.mat",...
%     "data_100604_9_Roberto.mat"];

gt_path = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Report_reduced_classification\Unified";
gt_peaks = ["ur_100526_5.mat",...
    "ur_100528_1.mat",...
    "ur_100528_4.mat",...
    "ur_100528_5.mat",...
    "ur_100528_6.mat",...
    "ur_100528_8.mat",...
    "ur_100601_1.mat",...
    "ur_100604_1.mat",...
    "ur_100604_9.mat"];

%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------


exp_name = regexprep(gta(i), '_Roberto.xlsx', '','ignorecase');
gts = readtable(fullfile(main_path,gta(i)));
gts = gts{:,2};
ats = load(fullfile(main_path,ata(i)));
at_cvector = ats.pclass;
load(fullfile(gt_path,gt_peaks(i)));

% Clean up 
for ii=1:1:length(gts)
    gts{ii} = regexprep(gts{ii}, 'unestable', 'Unknown','ignorecase');
    gts{ii} = regexprep(gts{ii}, 'non', 'no','ignorecase');
    gts{ii} = regexprep(gts{ii}, '\(.*?\)', '','ignorecase');    
    gts{ii} = regexprep(gts{ii}, ';', '+','ignorecase');
    gts{ii} = regexprep(gts{ii}, ':', '+','ignorecase');
    gts{ii} = regexprep(gts{ii}, 'oscillation[^s]','oscillations','ignorecase');
    gts{ii} = regexprep(gts{ii}, 'basal.+?no response','basal','ignorecase');
    gts{ii} = regexprep(gts{ii}, '\s*\+\s*', '+','ignorecase');
end
gts = lower(gts);

for ii=1:1:length(gts)    
    gts{ii} = regexprep(gts{ii},'basal activity\+.+','basal+response');
    gts{ii} = regexprep(gts{ii},'basal activity','basal');
    gts{ii} = regexprep(gts{ii},'peak.*','response');
    gts{ii} = regexprep(gts{ii},'oscillat.*','response');
    gts{ii} = regexprep(gts{ii},'plateau.*','response');
end

gt_cvector = zeros(length(gts),1);
for ii=1:1:length(gts)
    gt_cvector(ii) = find(clabels==gts{ii})-1;
end

cvector_u = unique([gt_cvector;at_cvector]);
for k=1:1:length(cvector_u)
    order_conf(k) = values(cdict,{cvector_u(k)});    
end 

gt_v = categorical(gt_cvector,cindex,clabels);
at_v = categorical(at_cvector,cindex,clabels);

C = confusionchart(gt_v,at_v);
C.ColumnSummary = 'column-normalized';
C.RowSummary = 'row-normalized';
ctit = regexprep(gta(i), '_','\\_','ignorecase');
C.Title = ctit;

hand_confusion = gca;
saveas(hand_confusion,fullfile(spath,"A_Confusion.png"));   

%% ROC
at_oh = onehotenc(cvector_u,at_cvector);
gt_oh = onehotenc(cvector_u,gt_cvector);
figure;
plotroc(gt_oh,at_oh)
hand_roc = gca;
saveas(hand_roc,fullfile(spath,"A_Roc.png")); 

%% profiles to plot
[ca_prof_smooth,~] = smoothProfiles(profiles_matrix,wfilter);
%% Differences
[difv] = vecdiff(at_cvector,gt_cvector);
figure;
hand_dif = gca;

disp(strcat("Looking for differences ..."));
disp(strcat("     ... 0 %"));
T_sample = 3;
t = (0:1:length([ats.cells_online(1).ca_profile])-1).*T_sample;
for k=1:1:length(difv)
    idx = difv(k);
    hold on
%     plot(t,ats.cells_online(idx).ca_profile,'LineWidth',1.5)
    plot(t,profiles_matrix(idx,:),'Color',[0.75 0.75 0.75],'LineWidth',1.5)    
    plot(t,ca_prof_smooth(idx,:),'b','LineWidth',1.5)
    ptitle=strcat("Profile ",num2str(idx),"  GT=",cdict(gt_cvector(idx)),"   Sys=",cdict(at_cvector(idx)));
    title(ptitle)
    xlabel('Time [s]')
    ylabel('Ratio')
    ylim([min([ats.cells_online.ca_profile]) max([ats.cells_online.ca_profile])])
    xlim([t(1) t(end)])
    xline(params(i,10)*T_sample,'--k','basal 0')
    xline(params(i,11)*T_sample,'--k','basal n')
    spk = [spks(idx).pks];
    if ~isempty(spk)
        [np, ~]=size(spk);
        for j=1:1:np
            plot(spk(j,1),spk(j,2),'or','MarkerSize',6)
            text(spk(j,1)-10,spk(j,2)+0.015,num2str(spk(j,2)),'Rotation',90,'Color','r')
        end
    end
    
    hold off   
    figname = strcat("B_diff",num2str(k,'%03d'),".png");
    saveas(hand_dif,fullfile(spath,figname));
    disp(strcat("     ... ",num2str(round(k*100/length(difv)))," %"));
    cla(hand_dif)
end
disp(strcat("Diffs Finished (:"));



%% Report generation
genReport(spath,exp_name)