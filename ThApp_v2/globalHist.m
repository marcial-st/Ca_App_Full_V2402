clear all;close all;clc

spath = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\ReportGen";
pathname = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Report_reduced_classification\Unified\090721";
filenames = ["100528_1.mat","100528_4.mat","100528_5.mat","100528_6.mat","100528_8.mat","100601_1.mat","100604_9.mat"];

prom_gt_global = [];
prom_all_global = [];

for k=1:length(filenames)
    if exist('profiles_smooth')==1
        clear profiles_smooth
    end
    if exist('spks')==1
        clear spks
    end
    disp(strcat("Reading ",filenames(k)))    
    disp("    progress...")
    load(fullfile(pathname,filenames(k)))
    [n_cells ~] = size(profiles_smooth);
    for i=1:1:n_cells
        %% Ground Truth
        [prom_gt,~,~,wact_gt,~] = getFeatures(spks(i).pks,profiles_smooth(i,:));
        prom_gt_global = [prom_gt_global;prom_gt];
        %% Plain Data
        [pks_all,pks_loc_all,pks_w_all,pks_prom_all] = findpeaks_wrapper(profiles_smooth(i,:));
%         [prom_all,~,~,wact_all,~] = getFeatures([pks_loc_all*Ts;pks_all]',profiles_smooth(i,:));
        prom_all_global = [prom_all_global;pks_prom_all];
        disp(strcat("    ",num2str(i/n_cells*100),"%"))
    end    
end
disp("Done!")

nbins = 200;
dist = 'normal';

figure;
h1 = histogram(prom_all_global,nbins)
h1.DisplayName = 'All peaks'
hold on
h2 = histogram(prom_gt_global,nbins)
h2.DisplayName = 'Ground Truth'
legend
title("Global Left Prominence")
ylim([0 500])
xlim([0 0.2])

figure;
subplot(1,2,1)
histfit(prom_all_global,nbins,dist)
subplot(1,2,2)
histfit(prom_gt_global,nbins,dist)