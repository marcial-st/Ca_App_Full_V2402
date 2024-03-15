function [class0,class1,class2,class3] = cmod_Basal_Response_ap2(profiles_matrix,t_start,visual_en,shft_l,shft_r,wfilter,ksec,t_basal_0,t_basal_n)
%%%%%%%%%%%%%%%%%%%%%%%%%% Vars & Params %%%%%%%%%%%%%%%%%%%%%%%%%%
% visual_en = 1;  % enable plots
[n_prof,n_frames] = size(profiles_matrix);
t = 1:1:n_frames;
%%%%%%%%%%%%%%%%%%%%%%%%%% Profile sets %%%%%%%%%%%%%%%%%%%%%%%%%%%

[ca_prof_smooth,noise_prof] = smoothProfiles(profiles_matrix,wfilter);

%%%%%%%%%%%%%%%%%%%%%%%% Classification %%%%%%%%%%%%%%%%%%%%%%%%%

% if visual_en;fprof = figure;end

% pk_h_th = 0.035;
pk_w_th = 0;
i_c0 = 1;
i_c1 = 1;
i_c2 = 1;
i_c3 = 1;
% shft = 10;
class0 =[];
class1 =[];
class2 =[];
class3 =[];

f = waitbar(0,'Cmod basalpeak progress...');
for i_prof = 1:n_prof  % Main for loop across profiles    
    
    pk_h_th = 6*nanstd(noise_prof(i_prof,t_basal_0:t_basal_n));
%     [pks,pks_loc,pks_w,pks_prom] = findpeaks(ca_prof_smooth(i_prof,shft_l:end-shft_r),'MinPeakProminence',pk_h_th,'WidthReference','halfprom');
    [pks,pks_prom] = islocalmax(ca_prof_smooth(i_prof,shft_l:end-shft_r),'MinSeparation',7,'MinProminence',0.0350);
    
    if isempty(pks)
        class0(i_c0,1) = i_prof;
        i_c0 = i_c0+1;
    else
        pks = pks(pks_w>=pk_w_th);
        if isempty(pks)
            class0(i_c0,1) = i_prof;
            i_c0 = i_c0+1;
        else
            pks_loc = pks_loc+shft_l-1;
            pks_loc = pks_loc(pks_w>=pk_w_th);
            pks_b = pks(pks_loc<t_start);
            pks_r = pks(pks_loc>=t_start);
            
            if (~isempty(pks_b)&&~isempty(pks_r))
                class3(i_c3,1) = i_prof;
                i_c3 = i_c3+1; 
            elseif(isempty(pks_b)&&~isempty(pks_r))
                class2(i_c2,1) = i_prof;
                i_c2 = i_c2+1;
            elseif(~isempty(pks_b)&&isempty(pks_r))
                class1(i_c1,1) = i_prof;
                i_c1 = i_c1+1;
            elseif(~isempty(pks_b)&&isempty(pks_r))
                class0(i_c0,1) = i_prof;
                i_c0 = i_c0+1;
            end                         
        end
    end         
%     if visual_en
%         figure(fprof);        
%         clf('reset');        
%         subplot(1,2,1)
%         hold on;
%         plot(t,ca_prof_smooth_plain(i_prof,:));         
%         plot(t(pks_loc),ca_prof_smooth_plain(i_prof,pks_loc),'*r')
%         hold off;
%         subplot(1,2,2)
%         plot(t,ca_prof_smooth(i_prof,:),'LineWidth',1); hold on;
%         plot(t,ca_prof_base(i_prof,:),'Color',[0.75 0.75 1]);
%         plot(t(pks_loc),ca_prof_smooth(i_prof,pks_loc),'*r')
%         ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
%         title(strcat("Prof ",num2str(i_prof)))
%         pause(0.2)          
%     end
    waitbar(i_prof/n_prof,f,'Cmod basalpeak progress...');
end
close(f)

if visual_en
    fclass0=figure;
    plot(t,ca_prof_smooth(class0,:)); hold on
    stem(t_start,max(max(ca_prof_smooth)),'r')
    title("Class no basal peak")
    ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])

    fclass1=figure;
    plot(t,ca_prof_smooth(class1,:)); hold on
    stem(t_start,max(max(ca_prof_smooth)),'r')
    title("Class with basal peaks")
    ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
end