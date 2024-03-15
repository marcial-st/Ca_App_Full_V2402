clear all;clc;close all;
%% MAIN IDX
i=2;
%% clean report path
% spath = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\ReportGen";
spath = "C:\Users\MST\Documents\ReportGen"
cleanSpath(spath)
%% Main 
cindex = [0,1,2,3,4];
clabels = ["unknown"; %-2
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
gt_path = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Peak_Analysis_Ground_Truth";

gt_a = ["100526_5_Ajelet.xlsx",...
    "100528_1_Ajelet.xlsx",...
    "100528_4_Ajelet.xlsx",...
    "100528_5_Ajelet.xlsx",...
    "100528_6_Ajelet.xlsx",...
    "100528_8_Ajelet.xlsx",...
    "100601_1_Ajelet.xlsx",...
    "100604_1_Ajelet.xlsx",...
    "100604_9_Ajelet.xlsx"];

gt_r = ["100526_5_Roberto.xlsx",...
    "100528_1_Roberto.xlsx",...
    "100528_4_Roberto_Activity.xlsx",...
    "100528_5_Roberto_Activity.xlsx",...
    "100528_6_Roberto_Activity.xlsx",...
    "100528_8_Roberto_Activity.xlsx",...
    "100601_1_Roberto_Activity.xlsx",...
    "100604_1_Roberto_Activity.xlsx",...
    "100604_9_Roberto_Activity.xlsx"];

% gt_peaks_a = ["data_100526_5_Ajelet.mat",...
%     "data_100528_1_Ajelet.mat",...
%     "data_100528_4_Ajelet.mat",...
%     "data_100528_5_Ajelet.mat",...
%     "data_100528_6_Ajelet.mat",...
%     "data_100528_8_Ajelet.mat",...
%     "data_100601_1_Ajelet.mat",...
%     "data_100604_1_Ajelet.mat",...
%     "data_100604_9_Ajelet.mat"];
% 
% gt_peaks_r = ["data_100526_5_Roberto.mat",...
%     "data_100528_1_Roberto.mat",...
%     "data_100528_4_Roberto.mat",...
%     "data_100528_5_Roberto.mat",...
%     "data_100528_6_Roberto.mat",...
%     "data_100528_8_Roberto.mat",...
%     "data_100601_1_Roberto.mat",...
%     "data_100604_1_Roberto.mat",...
%     "data_100604_9_Roberto.mat"];

gt_peaks_a = ["100526_5_Ajelet_Nuevo_test.mat",...
    "100528_1_Ajelet_Nuevo_test.mat",...
    "100528_4_Ajelet_Nuevo_test.mat",...
    "100528_5_Ajelet_Nuevo_test.mat",...
    "100528_6_Ajelet_Nuevo_test.mat",...
    "100528_8_Ajelet_Nuevo_test.mat",...
    "100601_1_Ajelet_Nuevo_test.mat",...
    "100604_1_Ajelet_Nuevo_test.mat",...
    "100604_9_Ajelet_Nuevo_test.mat"];

gt_peaks_r = ["100526_5_Roberto activity.mat",...
    "100528_1_Roberto activity.mat",...
    "100528_4_Roberto activity.mat",...
    "100528_5_Roberto activity.mat",...
    "100528_6_Roberto activity.mat",...
    "100528_8_Roberto activity.mat",...
    "100601_1_Roberto activity.mat",...
    "100604_1_Roberto activity.mat",...
    "100604_9_Roberto activity.mat"];

ata = ["ds_100526_5_r6.mat",...
    "ds_100528_1_r6.mat",...
    "ds_100528_4_r6.mat",...
    "ds_100528_5_r6.mat",...
    "ds_100528_6_r6.mat",...
    "ds_100528_8_r6.mat",...
    "ds_100601_1_r6.mat",...
    "ds_100604_4_r6.mat",...    
    "ds_100604_9_r6.mat"];

%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------


exp_name = regexprep(gt_r(i), '_Roberto.xlsx', '','ignorecase');

gts_a = readtable(fullfile(main_path,gt_a(i)));
gts_a = gts_a{:,2};

gts_r = readtable(fullfile(main_path,gt_r(i)));
gts_r = gts_r{:,2};

pks_a = load(fullfile(gt_path,gt_peaks_a(i)));
pks_r = load(fullfile(gt_path,gt_peaks_r(i)));

ats = load(fullfile(main_path,ata(i)));
at_cvector = ats.pclass;                                  


% Clean up 
for ii=1:1:length(gts_a)
    gts_a{ii} = regexprep(gts_a{ii}, 'unestable', 'Unknown','ignorecase');
    gts_a{ii} = regexprep(gts_a{ii}, 'non', 'no','ignorecase');
    gts_a{ii} = regexprep(gts_a{ii},  '\(.*?\)', '','ignorecase');    
    gts_a{ii} = regexprep(gts_a{ii}, ';', '+','ignorecase');
    gts_a{ii} = regexprep(gts_a{ii}, ':', '+','ignorecase');
    gts_a{ii} = regexprep(gts_a{ii}, 'oscillation[^s]','oscillations','ignorecase');
    gts_a{ii} = regexprep(gts_a{ii}, 'basal.+?no response','basal','ignorecase');
    gts_a{ii} = regexprep(gts_a{ii}, '\s*\+\s*', '+','ignorecase');
end
gts_a = lower(gts_a);

for ii=1:1:length(gts_a)    
    gts_a{ii} = regexprep(gts_a{ii},'basal activity\+.+','basal+response');
    gts_a{ii} = regexprep(gts_a{ii},'basal activity','basal');
    gts_a{ii} = regexprep(gts_a{ii},'peak.*','response');
    gts_a{ii} = regexprep(gts_a{ii},'oscillat.*','response');
    gts_a{ii} = regexprep(gts_a{ii},'plateau.*','response');
end

for ii=1:1:length(gts_r)
    gts_r{ii} = regexprep(gts_r{ii}, 'unestable', 'Unknown','ignorecase');
    gts_r{ii} = regexprep(gts_r{ii}, 'non', 'no','ignorecase');
    gts_r{ii} = regexprep(gts_r{ii},  '\(.*?\)', '','ignorecase');    
    gts_r{ii} = regexprep(gts_r{ii}, ';', '+','ignorecase');
    gts_r{ii} = regexprep(gts_r{ii}, ':', '+','ignorecase');
    gts_r{ii} = regexprep(gts_r{ii}, 'oscillation[^s]','oscillations','ignorecase');
    gts_r{ii} = regexprep(gts_r{ii}, 'basal.+?no response','basal','ignorecase');
    gts_r{ii} = regexprep(gts_r{ii}, '\s*\+\s*', '+','ignorecase');
end
gts_r = lower(gts_r);

for ii=1:1:length(gts_r)    
    gts_r{ii} = regexprep(gts_r{ii},'basal activity\+.+','basal+response');
    gts_r{ii} = regexprep(gts_r{ii},'basal activity','basal');
    gts_r{ii} = regexprep(gts_r{ii},'peak.*','response');
    gts_r{ii} = regexprep(gts_r{ii},'oscillat.*','response');
    gts_r{ii} = regexprep(gts_r{ii},'plateau.*','response');
end

%% -------------------------------------------
%%%                 Compare
% 
% if length(gts_a)==length(gts_r)
%     vdif = zeros(length(gts_a),1);
%     disp("Passed size check ")
%     disp("Compare check in progress ...")
%     
%     for k = 1:1:length(gts_a)        
%         if ~strcmp(gts_a{k},gts_r{k})
%             vdif(k) = 1;
%         end
%         disp(strcat("     ... ",num2str(round(k*100/length(gts_a)))," %"));
%     end
%     disp("     Compare finished (:")
% end
% 
% difidx = find(vdif==1);
% 
% for j=1:1:length(difidx)
%     disp(strcat(num2str(difidx(j))," Ajelet :",gts_a{difidx(j)}))
%     disp(strcat(num2str(difidx(j))," Roberto :",gts_r{difidx(j)}))
% end



%% Plots & Report
idx = 62;
Ts = 3;
wfilter = 11;
spath = "C:\Users\MST\Documents\ReportGen"
cleanSpath(spath)

profiles_matrix = pks_a.profiles_matrix;
[ca_prof_smooth,~] = smoothProfiles(profiles_matrix,wfilter);
t = [0:1:length(ats.cells_online(1).ca_profile)-1].*Ts;
figure;
hand_dif = gca;
disp(strcat("Ploting pks analysis ..."));
disp(strcat("     ... 0 %"));
i_lim = length(ats.cells_online);
for idx=1:1:i_lim
    
    hold on
    plot(t,profiles_matrix(idx,:),'Color',[0.65 0.65 0.65],'LineWidth',1.5,'DisplayName','Original')
    plot(t,ca_prof_smooth(idx,:),'b','LineWidth',1.5,'DisplayName','Smooth')
    ptitle=strcat("Profile ",num2str(idx));
    title(ptitle)
    xlabel('Time [s]')
    ylabel('Ratio')
    ylim([min([ats.cells_online.ca_profile]) max([ats.cells_online.ca_profile])])
    xlim([t(1) t(end)])
    spk_a = [pks_a.spks(idx).pks];
    if ~isempty(spk_a)
      [npa, ~]=size(spk_a);
      for j=1:1:npa
%          plot(spk_a(j,1),spk_a(j,2),'m+','MarkerSize',10,'DisplayName','A','LineWidth',1.5)
         plot(spk_a(j,1),spk_a(j,2),'m+','MarkerSize',10,'LineWidth',1.5)
         text(spk_a(j,1)-10,spk_a(j,2)+0.015,num2str(spk_a(j,2)),'Rotation',90,'Color','m')
      end
      text(20,max([ats.cells_online.ca_profile])-0.015,'+ A','Color','m','FontSize',11)
    end
    spk_r = [pks_r.spks(idx).pks];
    if ~isempty(spk_r)
      [npr, ~]=size(spk_r);
      for j=1:1:npr
%         plot(spk_r(j,1),spk_r(j,2),'ok','MarkerSize',7,'DisplayName','R','LineWidth',1)
        plot(spk_r(j,1),spk_r(j,2),'ok','MarkerSize',7,'LineWidth',1)
        text(spk_r(j,1)+10,spk_r(j,2)+0.015,num2str(spk_r(j,2)),'Rotation',90,'Color','k')
      end
      text(20,max([ats.cells_online.ca_profile])-0.035,'o R','Color','k','FontSize',11)
    end    
    legend('Original','Smooth')
%     set(gca,'Color',[0.7 0.7 0.7])
    hold off
%     pause(0.2)
    
    figname = strcat("B_diff",num2str(idx,'%03d'),".png");
    saveas(hand_dif,fullfile(spath,figname));

    if idx ~= i_lim
        cla(hand_dif)
    end
    
    disp(strcat("     ... ",num2str(round(idx*100/length(ats.cells_online)))," %"));
end
disp(" Plot finished (:")
%% Report generation
genReport(spath,exp_name)