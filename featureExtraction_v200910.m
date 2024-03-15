% function [CELLSf] = featureExtraction_v200910(profiles_matrix)
clear all;clc
aux = load('online_profiles.mat');
profiles_matrix = aux.online_profiles';

zcIdx = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %function to return zero indices
[n_prof,n_frames] = size(profiles_matrix);


ca_prof_smooth = zeros(n_prof,n_frames);   %% Smooth profiles with moving average filter
for i_smooth = 1:1:n_prof
    ca_prof_smooth(i_smooth,:) = smooth(profiles_matrix(i_smooth,:),5);
end

ca_prof_diff = round(diff(ca_prof_smooth')',6);  %% Ca Profs differences matrix
max_diff = max(max(ca_prof_diff));
[~,t_start] = find(ca_prof_diff == max_diff);
t_start = t_start-4;

n_segments = 3;
n_segvector = floor(linspace(t_start,n_frames,n_segments+1));

struct_size = n_prof;  %number of elements
CELLSf(struct_size) = struct();
thMainPeakProm = 0.04;
thMainPeakH = 0.1;
% f = waitbar(0,'Please wait...');
figure;
% for i_prof=1:1:n_prof
i_prof = 120;
    flag_mainpk = 0;
    ca_prof_i = ca_prof_smooth(i_prof,:);
    CELLSf(i_prof).Basal = nanmean(ca_prof_i(1:t_start));
    
    ca_prof_unbias = ca_prof_i-CELLSf(i_prof).Basal;
    [pks,pks_l] = findpeaks(ca_prof_unbias(t_start:end),'MinPeakProminence',0.01); %3*std(ca_prof_unbias)
    CELLSf(i_prof).nPeaks = length(pks);
    
    plot(ca_prof_unbias); hold on;
    plot(pks_l+t_start,pks,'m+')
    title(strcat("Profile:",num2str(i_prof)))
    
    for i_seg = 1:length(n_segvector)-1
        ca_prof_seg = ca_prof_unbias(n_segvector(i_seg)+1:n_segvector(i_seg+1));
        t=n_segvector(i_seg)+[1:1:length(ca_prof_seg)];
        [pks,pks_locs,pks_width,pks_prom] = findpeaks(ca_prof_seg,'MinPeakProminence',0.02,'Annotate','extents','WidthReference','halfheight'); %3*std(ca_prof_unbias)        
        [vly,vly_locs,~,vly_prom] = findpeaks(-ca_prof_seg,'MinPeakProminence',0.005);
        if ~isempty(pks)
            CELLSf(i_prof).(strcat('pks_val_',num2str(i_seg))) = mean(pks);
            CELLSf(i_prof).(strcat('pks_num_',num2str(i_seg))) = length(pks);
            CELLSf(i_prof).(strcat('pks_width_',num2str(i_seg))) = mean(pks_width);
        else
            CELLSf(i_prof).(strcat('pks_val_',num2str(i_seg))) = 0;
            CELLSf(i_prof).(strcat('pks_num_',num2str(i_seg))) = 0;
            CELLSf(i_prof).(strcat('pks_width_',num2str(i_seg))) = 0;
        end
                 
        plot(pks_locs+n_segvector(i_seg),pks,'ob')
        plot(n_segvector,zeros(length(n_segvector)),'*k')
        
        if (~isempty(pks) && (flag_mainpk == 0) && (i_seg == 1))
            pk_relevant = find(pks_prom>thMainPeakProm);
            if isempty(pk_relevant)
                pk_relevant = find(pks>thMainPeakH);
            end
            if ~(isempty(pk_relevant))
                CELLSf(i_prof).MainPeakVal = ca_prof_seg(pks_locs(pk_relevant(1)))+CELLSf(i_prof).Basal; %% Feature Main Peak Value
                CELLSf(i_prof).MainPeakHeight = ca_prof_seg(pks_locs(pk_relevant(1))); %% Feature Main Peak Prominence    
                CELLSf(i_prof).MainPeakWidth = pks_width(pk_relevant(1));                   %% Feature Main Peak Width
%                 CELLSf(i_prof).PeakTime = pks_locs(pk_relevant(1))+t_start;
                CELLSf(i_prof).PeakTime = pks_locs(pk_relevant(1))+n_segvector(i_seg);
                
                plot(vly_locs+n_segvector(i_seg),ca_prof_seg(vly_locs),'or')
                vly_locs = vly_locs(vly_locs>pks_locs(pk_relevant(1)));
                if isempty(vly_locs)
                    ca_down_seg = ca_prof_seg(pks_locs(pk_relevant(1)):end);
                else                    
                    ca_down_seg = ca_prof_seg(pks_locs(pk_relevant(1)):vly_locs(1));
                end
                
                mpk90_e = 0.9*CELLSf(i_prof).MainPeakHeight;
                mpk60_e = 0.6*CELLSf(i_prof).MainPeakHeight;                
                
                [mpk90,fTimePk90] = min(abs(ca_down_seg-mpk90_e));
                [mpk60,fTimePk60] = min(abs(ca_down_seg-mpk60_e));
                
                if abs(mpk90_e-ca_down_seg(fTimePk90))<=0.02
                    CELLSf(i_prof).MainPk90 = ca_down_seg(fTimePk90);
                    CELLSf(i_prof).Main_T_Pk90 = fTimePk90(1)+CELLSf(i_prof).PeakTime;
                    plot(CELLSf(i_prof).Main_T_Pk90,CELLSf(i_prof).MainPk90,'r*')
                else
                    CELLSf(i_prof).MainPk90 = 0;
                    CELLSf(i_prof).Main_T_Pk90 = 0;
                end
                
                if abs(mpk60_e-ca_down_seg(fTimePk60))<=0.02
                    CELLSf(i_prof).MainPk60 = ca_down_seg(fTimePk60);
                    CELLSf(i_prof).Main_T_Pk60 = fTimePk60(1)+CELLSf(i_prof).PeakTime;
                    plot(CELLSf(i_prof).Main_T_Pk60,CELLSf(i_prof).MainPk60,'r*')
                else
                    CELLSf(i_prof).MainPk60 = 0;
                    CELLSf(i_prof).Main_T_Pk60 = 0;
                end
                
                plot(CELLSf(i_prof).PeakTime,CELLSf(i_prof).MainPeakHeight,'r*')                                
                
                flag_mainpk = 1;
            else
                CELLSf(i_prof).MainPeakVal = 0;
                CELLSf(i_prof).MainPeakHeight = 0;
                CELLSf(i_prof).MainPeakWidth = 0;
                CELLSf(i_prof).PeakTime = 0;
              
                CELLSf(i_prof).MainPk90 = 0;
                CELLSf(i_prof).MainPk60 = 0;                
                
                CELLSf(i_prof).Main_T_Pk90 = 0;
                CELLSf(i_prof).Main_T_Pk60 = 0;                
            end
        else
            CELLSf(i_prof).MainPeakVal = 0;
            CELLSf(i_prof).MainPeakHeight = 0;
            CELLSf(i_prof).MainPeakWidth = 0;
            CELLSf(i_prof).PeakTime = 0;
              
            CELLSf(i_prof).MainPk90 = 0;
            CELLSf(i_prof).MainPk60 = 0;            
                
            CELLSf(i_prof).Main_T_Pk90 = 0;
            CELLSf(i_prof).Main_T_Pk60 = 0;                       
        end
        
%         if i_seg ~= 1
%            m=t'\ca_prof_seg'; 
%            y1 = m*t; 
%            plot(t,y1,'m')
           p = polyfit(t,ca_prof_seg,2);
           y2 = polyval(p,t);
           plot(t,y2,'g')
           CELLSf(i_prof).(strcat('fitc1_',num2str(i_seg)))=p(1);
           CELLSf(i_prof).(strcat('fitc2_',num2str(i_seg)))=p(2);
           CELLSf(i_prof).(strcat('fitc3_',num2str(i_seg)))=p(3);
%         end
        
        CELLSf(i_prof).(strcat('auc_',num2str(i_seg)))= trapz(ca_prof_unbias(n_segvector(i_seg)+1:n_segvector(i_seg+1)));
        CELLSf(i_prof).(strcat('kurt_',num2str(i_seg)))= kurtosis(ca_prof_unbias(n_segvector(i_seg)+1:n_segvector(i_seg+1)));
        CELLSf(i_prof).(strcat('entr_',num2str(i_seg)))= entropy(ca_prof_unbias(n_segvector(i_seg)+1:n_segvector(i_seg+1))/max(abs(ca_prof_unbias(n_segvector(i_seg)+1:n_segvector(i_seg+1)))));
        CELLSf(i_prof).(strcat('zcross_',num2str(i_seg)))= length(zcIdx(ca_prof_unbias(n_segvector(i_seg)+1:n_segvector(i_seg+1))));
    end
    hold off;
%     waitbar(i_prof/n_prof,f,'Extracting Features');
%     pause(0.1)        
% end
% close(f)