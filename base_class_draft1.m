function [class_1_peak,class_0_peak] = base_class_draft1(profiles_matrix)
% clear all;close all;clc
% load prof_mat_1.mat
tic
visual_en = 0;
eps_val = 0.01;
zcIdx = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %function to return zero indices
[n_prof,n_frames] = size(profiles_matrix);
t = 1:1:n_frames;

ca_prof_smooth = zeros(n_prof,n_frames);   %% Smooth profiles with moving average filter
% if visual_en == 1
%     fprof = figure;
% end
for i_smooth = 1:1:n_prof
    s_prof = smooth(profiles_matrix(i_smooth,:),5)';
    p_min = islocalmin(s_prof);
    p_min(1)=1;p_min(end)=1; 
    ca_prof_base(i_smooth,:) = interp1(t(p_min),s_prof(p_min),t,'linear');        
    ca_prof_smooth(i_smooth,:) = s_prof;
    ca_prof_smooth_plain(i_smooth,:) = s_prof-ca_prof_base(i_smooth,:);
%     if visual_en == 1        
%         figure(fprof)
%         plot(t,ca_prof_smooth(i_smooth,:),t,ca_prof_smooth_plain(i_smooth,:)+mean(ca_prof_smooth(i_smooth,1:10)),'r',t,ca_prof_base(i_smooth,:),'g');
%         ylim([0,1.5])
%         xlim([t(1),t(end)])
%         title(strcat("Profile: ",num2str(i_smooth),"/",num2str(n_prof)))
%         pause(0.2)
%     end    
end


ca_prof_diff = round(diff(ca_prof_smooth_plain')',6);  %% Ca Profs differences matrix
max_diff = max(max(ca_prof_diff));
[~,t_start] = find(ca_prof_diff == max_diff);
t_start = t_start-8;
t_start = 60;
% if visual_en == 1
    fstart = figure;
    figure(fstart)
    plot(t,profiles_matrix')
    hold on;
    stem(t_start,1.5,'LineWidth',1.5,'color','m')
    xlim([t(1),t(end)])
    title("Base classification")
    hold off;
% end

pk_th = 10*mean(nanstd(ca_prof_smooth_plain(:,1:t_start),0,2));

struct_size = n_prof;  %number of elements
CELLSf(struct_size) = struct();
class_count = 1;
n_segments = 3;
n_segvector = floor(linspace(t_start,n_frames,n_segments+1));


for i_prof=1:1:n_prof
% i_prof = 70;
    flag_mainpk = 0;
    
    CELLSf(i_prof).Basal = nanmean(ca_prof_smooth(i_prof,1:t_start));
    
    ca_prof_unbias = ca_prof_smooth(i_prof,:)-CELLSf(i_prof).Basal;    
    ca_prof_plain_unbias = ca_prof_smooth_plain(i_prof,:)-nanmean(ca_prof_smooth_plain(i_prof,1:t_start));
    
%     if visual_en == 1
%         figure(main_fig)
%         plot(ca_prof_unbias,'b','LineWidth',1,'DisplayName','UnbCaProf'); legend; hold on;
%         stem(t_start,max(ca_prof_unbias)+0.05,'Color',[0.8,0.8,0.8],'DisplayName','tstart');legend
%         plot(n_segvector,zeros(length(n_segvector)),'+k');
%         title(strcat("Profile:",num2str(i_prof)))
%     end
    
    [pks,pks_l] = findpeaks(ca_prof_plain_unbias(1:100),'MinPeakHeight',pk_th);
    [~,vly_l] = findpeaks(-ca_prof_plain_unbias(1:100)+max(ca_prof_plain_unbias),'MinPeakHeight',pk_th);
    
    pks_1_basal = pks_l(pks_l<t_start);
    pks_l_sec = pks_l(pks_l>=t_start);
    
    t_main_pk =[];
    t_main_end =[];
    
    
    if ~isempty(pks)
        if visual_en == 1  
            main_fig = figure;
            figure(main_fig)
            plot(ca_prof_unbias,'b','LineWidth',1,'DisplayName','UnbCaProf'); legend; hold on;
            stem(t_start,max(ca_prof_unbias)+0.05,'Color',[0.8,0.8,0.8],'DisplayName','tstart');legend
            plot(n_segvector,zeros(length(n_segvector)),'+k');
            title(strcat("Profile:",num2str(i_prof)))
        end
        
        CELLSf(i_prof).nPeaksBasal = length(pks_1_basal);     %% Basal peaks        
        if ~isempty(pks_1_basal)
            if visual_en == 1  
                figure(main_fig)
                plot(t(pks_1_basal),ca_prof_unbias(pks_1_basal),'r*','MarkerSize',8,'DisplayName','PeakBasal')
            end
        end
        
        if ((length(pks_l_sec)>=1)&&(~isempty(vly_l)))
            t_main_pk= pks_l_sec(1);
            aux_locs = vly_l(vly_l>t_main_pk);
            class_1_peak(class_count) = i_prof;
            class_count = class_count+1;
            if ~isempty(aux_locs)
                t_main_end = aux_locs(1);
                if visual_en == 1;stem(t_main_end,max(ca_prof_unbias)+0.05,'Color',[1,0.8,0.8]);end
            end         
        end
                        
        if length(pks_l_sec)>1                                %% Secondary peaks
            CELLSf(i_prof).nPeaksSec = length(pks_l_sec)-1;
            if visual_en ==1
                figure(main_fig);
                plot(pks_l_sec(1),ca_prof_unbias(pks_l_sec(1)),'r*','MarkerSize',9,'DisplayName','PeakMain');
                plot(pks_l_sec(2:end),ca_prof_unbias(pks_l_sec(2:end)),'m+','MarkerSize',7,'DisplayName','PeaksSec');                
                legend;                
            end
        else
            CELLSf(i_prof).nPeaksSec = length(pks_l_sec);
        end 
    else
        CELLSf(i_prof).nPeaksBasal = 0;
        CELLSf(i_prof).nPeaksSec = 0;
    end
                               

    hold off;
    waitbar(i_prof/n_prof,f,'Extracting Features');
%     pause(0.3)        
end
class_0_peak = 1:1:n_prof;
class_0_peak(class_1_peak)=[];

% figure;plot(t,ca_prof_smooth_plain(class_1_peak,:));title("Class with peaks");ylim([-0.1,0.5])
% figure;plot(t,ca_prof_smooth_plain(class_0_peak,:));title("Class without peaks");ylim([-0.1,0.5])

figure;plot(t,ca_prof_smooth(class_1_peak,:));title("Class with peaks");ylim([0.5,1.5])
figure;plot(t,ca_prof_smooth(class_0_peak,:));title("Class without peaks");ylim([0.5,1.5])

close(f)
toc