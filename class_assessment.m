clear;clc;
cindex = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
clabels = ["Unknown"; %0
    "No response"; %1
    "Peak"; %2
    "Peak+plateau"; %3
    "Peak+oscillations";%4
    "Peak+plateau+oscillations";%5    
    "Oscillations";%6   
    "Oscillations+plateau";%7
    "Basal activity+Peak";%8
    "Basal activity+Peak+plateau";%9
    "Basal activity+Peak+oscillations";%10
    "Basal activity+Peak+plateau+oscillations";%11
    "Basal activity+Oscillations"; %12
    "Basal activity+Oscillations+plateau"; %13
    "Basal activity";%14
    "Basal activity+plateu";%15
    ];
clabels = lower(clabels);
cdict = containers.Map(cindex,clabels);
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
main_path = "C:\Users\MST\OneDrive\Perfiles\PeakGT";
gta = ["100526_5_Roberto.xlsx",...
    "100528_1_Roberto.xlsx",...
    "100528_4_Roberto.xlsx",...
    "100528_5_Roberto.xlsx",...
    "100528_6_Roberto.xlsx",...
    "100528_8_Roberto.xlsx",...
    "100601_1_Roberto.xlsx",...
    "100604_1_Roberto.xlsx",...
    "100604_9_Roberto.xlsx"];
   
ata = ["ds_100526_5_2.mat",...
    "ds_100528_1_2.mat",...
    "ds_100528_4_2.mat",...
    "ds_100528_5_2.mat",...
    "ds_100528_6_2.mat",...
    "ds_100528_8_2.mat",...
    "ds_100601_1_2.mat",...
    "ds_100604_4_2.mat",...    
    "ds_100604_9_2.mat"];
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
i=5;

gts = readtable(fullfile(main_path,gta(i)));
gts = gts{:,2};
ats = load(fullfile(main_path,ata(i)));
at_cvector = ats.pclass;




% Clean up 
for ii=1:1:length(gts)
    gts{ii} = regexprep(gts{ii}, 'unestable', 'Unknown','ignorecase');
    gts{ii} = regexprep(gts{ii}, 'non', 'no','ignorecase');
    gts{ii} = regexprep(gts{ii}, '\(.*?\)', '','ignorecase');    
    gts{ii} = regexprep(gts{ii}, ';', '+','ignorecase');
    gts{ii} = regexprep(gts{ii}, ':', '+','ignorecase');
    gts{ii} = regexprep(gts{ii}, 'oscillation[^s]','oscillations','ignorecase');
    gts{ii} = regexprep(gts{ii}, 'basal.+?no response','basal activity','ignorecase');
    gts{ii} = regexprep(gts{ii}, '\s*\+\s*', '+','ignorecase');
end
gts = lower(gts);

gt_cvector = zeros(length(gts),1);
for ii=1:1:length(gts)
    gt_cvector(ii) = find(clabels==gts{ii})-1;
end

cvector_u = unique([gt_cvector;at_cvector]);
for k=1:1:length(cvector_u)
    order_conf(k) = values(cdict,{cvector_u(k)});    
end 

C = confusionchart(gt_cvector,at_cvector);
C.ColumnSummary = 'column-normalized';
C.RowSummary = 'row-normalized';
ctit = regexprep(gta(i), '_','\\_','ignorecase');
C.Title = ctit;

%% ROC
at_oh = onehotenc(cvector_u,at_cvector);
gt_oh = onehotenc(cvector_u,gt_cvector);
figure;
plotroc(gt_oh,at_oh)


%% Differences
[difv] = vecdiff(at_cvector,gt_cvector);

idx = difv(1);

figure;
plot(ats.cells_online(idx).ca_profile)
ptitle=strcat("Profile ",num2str(idx),"  GT=",cdict(gt_cvector(idx)),"   Sys=",cdict(at_cvector(idx)));
title(ptitle)