%% Script to identify and edit GT datasets 
clear all;close all; 

spath = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\ReportGen";
pathname = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Report_reduced_classification\Unified\111721_full_basal";
filenames = ["100528_1.mat","100528_4.mat","100528_5.mat","100528_6.mat","100528_8.mat","100601_1.mat","100604_9.mat"];
tstims = [90,150,126,192,189,189,189];
epsilon = 10;
disp(strcat("Removing redundancy within epsilon = ",num2str(epsilon)))
Ts = 3;
for kf = 1:1:length(filenames)
    filename = filenames(kf);
    disp(strcat("Fixing Experiment: ",filename))
    data = load(fullfile(pathname,filename));
    [np,nf] = size(data.profiles_smooth);
    figure;
    time = [0:1:nf-1].*Ts;
    profiles = data.profiles_smooth;
    profiles = profiles(60:90,:);
    basals = mean(profiles(:,1:30),2);
    profiles_aligned = profiles - basals;
    plot(time,profiles_aligned);
    title(filename,"FontSize",20)
    xlim([time(1) time(end)])    
    ylabel("Ratio")
    xlabel("time [s]")
%     line_label = strcat("t_{stim}=",tstims(kf));
    xline(tstims(kf),'--',{strcat('t_{stim}=',num2str(tstims(kf)))},'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','LabelOrientation','horizontal','LineWidth',2)
%     if isfield(data,'spks')
%         f = waitbar(0,strcat("Exp. ",filename," progress"));
%         pause(0.5)
%         for i=1:1:length(data.spks)
%             waitbar(i/length(data.spks),f,strcat("Exp. ",filename," progress = ",num2str(i),"/",num2str(length(data.spks))));
%             aux = [data.spks(i).pks];
%             if ~isempty(aux)
%                 [pk_max_loc_f,pks_all] = islocalmax(data.profiles_smooth(i,:));
%                 locs = find(pk_max_loc_f);
%                 [max_j,~]=size(aux);
%                 for j=1:1:length(max_j)
%                     loc_error = abs(locs-aux(j,1));
%                     [~,idx] = min(loc_error);
%                     if min(loc_error)>0
%                         disp(strcat("    Prof.",num2str(i)," Original GT @",num2str(aux(j,1)),"  updated to @",num2str(locs(idx))))
%                     end
%                     aux(j,1) = locs(idx);
%                 end
%                 data.spks(i).pks = aux;
%             end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %               Block to remove GT redundancy            
% %             if ~isempty(aux)             
% %                 idx = 1;
% %                 dlength = length(aux);
% %                 while idx<dlength
% %                     aux2 = abs(aux(:,1)-aux(idx,1));
% %                     idx2 = find(aux2<=epsilon);
% %                     if length(idx2)>1
% %                         disp(strcat("Warning: Duplicated GT in profile=",num2str(i)," location="));                        
% %                         aux(idx2,:)
% %                         aux(idx2(2:end),:)=[];                        
% %                     end
% %                     new_pks(idx,:) = aux(idx2(1),:);
% %                     idx = idx+1;
% %                     dlength = length(aux);
% %                 end
% %                 data.spks(i).pks = aux;
% %             end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end
%         save(fullfile(pathname,filename),'-struct','data')
%         disp("...saved")
%         close(f)
%     else
%         disp(strcat("Warning! Experiment ",filename," does not contain spks structure"))
%     end
    disp("...done")
end