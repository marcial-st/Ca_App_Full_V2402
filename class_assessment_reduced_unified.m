clear all;clc;close all;
%% MAIN IDX
i=9;
%% clean report path
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
gt_path = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Report_reduced_classification\Unified";

gt_peaks_u = ["ur_100526_5.mat",...
    "ur_100528_1.mat",...
    "ur_100528_4.mat",...
    "ur_100528_5.mat",...
    "ur_100528_6.mat",...
    "ur_100528_8.mat",...
    "ur_100601_1.mat",...
    "ur_100604_9.mat"];
% "ur_100604_1.mat",...

ata = ["ds_100526_5_r6.mat",...
    "ds_100528_1_r6.mat",...
    "ds_100528_4_r6.mat",...
    "ds_100528_5_r6.mat",...
    "ds_100528_6_r6.mat",...
    "ds_100528_8_r6.mat",...
    "ds_100601_1_r6.mat",...    
    "ds_100604_9_r6.mat"];
%   "ds_100604_4_r6.mat",...

%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
wb = waitbar(0,'Overall progress, please graba a coffee...');
prom_gl = [];
prom_gr = [];
prom_gmax = [];
wact_g = [];
dstim_g =[];
class_g = [];
prom_gln = [];
wact_gn = [];
dstim_gn =[];
class_gn = [];

for i=1:1:5

exp_name = regexprep(gt_peaks_u(i), '.mat', '','ignorecase');
exp_name = regexprep(exp_name, 'ur_', '','ignorecase');

pks_u = load(fullfile(gt_path,gt_peaks_u(i)));

ats = load(fullfile(main_path,ata(i)));
at_cvector = ats.pclass;                                  

%% Plots & Report
% idx = 62;
Ts = 3;
wfilter = 11;
spath = "C:\Users\MST\Documents\ReportGen"
cleanSpath(spath)

profiles_matrix = pks_u.profiles_matrix;
[ca_prof_smooth,~] = smoothProfiles(profiles_matrix,wfilter);
t = [0:1:length(ats.cells_online(1).ca_profile)-1].*Ts;
hand_dif = figure;

% hand_dif = gca;
disp(strcat("Ploting pks analysis ..."));
disp(strcat("     ... 0 %"));
i_lim = length(ats.cells_online);
% prom_gl = [];
% prom_gr = [];
% prom_gmax = [];
% wact_g = [];
for idx=1:1:i_lim
    
    subplot(2,1,1)
    hold on
    plot(t,profiles_matrix(idx,:),'Color',[0.65 0.65 0.65],'LineWidth',1.5,'DisplayName','Original')
    plot(t,ca_prof_smooth(idx,:),'b','LineWidth',1.5,'DisplayName','Smooth')
    ptitle=strcat("Profile ",num2str(idx));
    title(ptitle)
    xlabel('Time [s]')
    ylabel('Ratio')
    ylim([min([ats.cells_online.ca_profile]) max([ats.cells_online.ca_profile])])
    xlim([t(1) t(end)])
    spk_u = [pks_u.spks(idx).pks];
    if ~isempty(spk_u)
        %% sort spk_u
        [~,ik] = sort(spk_u(:,1));
        spk_u = spk_u(ik(:,1),:);
        %%        
      [npa, ~]=size(spk_u);
%       [prom_l,prom_r,prom_max,wact,spk_u,dstim] = getProminence(spk_u,ca_prof_smooth,tstim)
      [prom_l,prom_r,prom_max,wact,spk_u,dstim,class] = getProminence(spk_u,ca_prof_smooth(idx,:),pks_u.params.tstim,0);
      prom_gl = [prom_gl;prom_l];
      prom_gr = [prom_gr;prom_r];
      prom_gmax = [prom_gmax;prom_max];
      wact_g = [wact_g;wact];
      dstim_g = [dstim_g;dstim];
      class_g = [class_g;class];
      
      spk_max = 3*find(islocalmax(ca_prof_smooth(idx,:),'MinSeparation',7,'MinProminence',0.0075))';
      spk_max = cleanmax(spk_max,spk_u(:,1));
      [prom_n,~,~,wn,~,dstimn,classn] = getProminence(spk_max,ca_prof_smooth(idx,:),pks_u.params.tstim,1);
      prom_gln = [prom_gln;prom_n];
      wact_gn = [wact_gn;wn];
      dstim_gn = [dstim_gn;dstimn];
      class_gn = [class_gn;classn];
      for j=1:1:npa
%          plot(spk_u(j,1),spk_u(j,2),'m+','MarkerSize',10,'DisplayName','A','LineWidth',1.5)         
         plot(spk_u(j,3),spk_u(j,4),'m+','MarkerSize',10,'LineWidth',1.5)
         text(spk_u(j,3)-10,spk_u(j,4)+0.015,num2str(prom_l(j)),'Rotation',90,'Color','m')
      end
      text(20,max([ats.cells_online.ca_profile])-0.015,'+ A','Color','m','FontSize',11)
    end
%     set(gca,'Color',[0.7 0.7 0.7])

    legend('Original','Smooth')
    hold off

    subplot(2,1,2)
    if ~isempty(spk_u)
        histogram(prom_l,'NumBins',20)
        xlabel('Prom Value')
        ylabel('Freq')
    end
%     pause(0.2)
    figname = strcat("B_diff",num2str(idx,'%03d'),".png");
    saveas(hand_dif,fullfile(spath,figname));

    if idx ~= i_lim
        clf(hand_dif)
    end
    
    disp(strcat("     ... ",num2str(round(idx*100/length(ats.cells_online)))," %"));
end

hedges = 0:0.02:0.4;

ghist = figure;
subplot(3,1,1)
hpl = histogram(prom_gl,'BinEdges',hedges)
title('Left prominence')
xlabel('Prom Value')
ylabel('Freq')
x = hpl.BinEdges ;
Prom_L = hpl.Values ;
text(x(1:end-1),Prom_L,num2str(Prom_L'),'vert','bottom','horiz','center');
subplot(3,1,2)
hpr = histogram(prom_gr,'BinEdges',hedges)
title('Rigth prominence')
xlabel('Prom Value')
ylabel('Freq')
x = hpr.BinEdges ;
Prom_R = hpr.Values ;
text(x(1:end-1),Prom_R,num2str(Prom_R'),'vert','bottom','horiz','center');
subplot(3,1,3)
hpm = histogram(prom_gmax,'BinEdges',hedges)
title('Max prominence')
xlabel('Prom Value')
ylabel('Freq')
x = hpm.BinEdges ;
Prom_max = hpm.Values ;
text(x(1:end-1),Prom_max,num2str(Prom_max'),'vert','bottom','horiz','center');
saveas(ghist,fullfile(spath,"B_bglobalhist.png"));

Ranges=[strcat(num2str(x(1:end-1)','%0.02f'),"-",num2str(x(2:end)','%0.02f'))];
mtable = table(Ranges,Prom_L',Prom_R',Prom_max');
mtable.Properties.VariableNames = {'Ranges','Prom_L','Prom_R','Prom_max'}

data.min_prom_gl = min(nonzeros(prom_gl));
data.max_prom_gl = max(prom_gl);
data.min_prom_gr = min(nonzeros(prom_gr));
data.max_prom_gr = max(prom_gr);
% data.min_prom_gmax = min(nonzeros(prom_gmax));
% data.max_prom_gmax = max(prom_gmax);

fwidth=figure;
wedges = 0:5:40;
hwact = histogram(wact_g,'BinEdges',wedges,'FaceColor','g')
title('Activity Width')
xlabel('Frames')
ylabel('Freq')
x = hwact.BinEdges ;
y = hwact.Values ;
text(x(1:end-1),y,num2str(y'),'vert','bottom','horiz','center');
saveas(fwidth,fullfile(spath,"B_bglobalw.png"));

waitbar(i/8,wb,'Overall progress, please graba a coffee...');
end
% close(wb)
disp(" Plot finished (:")

%% Report generation
%  genReport(spath,exp_name,mtable,data)
%%
zidx = find(wact_g==0);
prom_gl(zidx)=[];
wact_g(zidx)=[];
dstim_g(zidx)=[];
class_g(zidx)=[];
class_g = char(class_g);

zidx = find(wact_gn==0);
prom_gln(zidx)=[];
wact_gn(zidx)=[];
dstim_gn(zidx)=[];
class_gn(zidx)=[];
class_gn = char(class_gn);

vaux = 1:1:length(prom_gln);
vaux = randsample(vaux,length(prom_gl));

feat_act = table(prom_gl,wact_g,dstim_g,class_g,'VariableNames',{'prominence','width','distance','class'});
feat_nact = table(prom_gln,wact_gn,dstim_gn,class_gn,'VariableNames',{'prominence','width','distance','class'});
feat_nactr = feat_nact(vaux,:);
features = [feat_act;feat_nactr];

%% Region1: Basal
idx_r1 = find(feat_act.distance<=0);
feat_act_r1 = feat_act(idx_r1,:);

idx_r1n = find(feat_nact.distance<=0);
vaux = randsample(idx_r1n,length(idx_r1));
feat_nact_r1 = feat_nact(vaux,:);
features_r1 = [feat_act_r1;feat_nact_r1];

%% Region2: No Basal
idx_r2 = find(feat_act.distance>0);
feat_act_r2 = feat_act(idx_r2,:);

idx_r2n = find(feat_nact.distance>0);
vaux = randsample(idx_r2n,length(idx_r2));
feat_nact_r2 = feat_nact(vaux,:);
features_r2 = [feat_act_r2;feat_nact_r2];
disp('Done!')