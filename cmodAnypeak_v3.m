function [class1,class0] = cmodAnypeak_v3(profiles_matrix,visual_en,shft_l,shft_r,wfilter,ksec,t_start,t_basal_0,t_basal_n)
%%%%%%%%%%%%%%%%%%%%%%%%%% Vars & Params %%%%%%%%%%%%%%%%%%%%%%%%%%
% visual_en = 1;  % enable plots
[n_prof,n_frames] = size(profiles_matrix);
t = 1:1:n_frames;
%%%%%%%%%%%%%%%%%%%%%%%%%% Profile sets %%%%%%%%%%%%%%%%%%%%%%%%%%%

[ca_prof_smooth,noise_prof] = smoothProfiles(profiles_matrix,wfilter);

%%%%%%%%%%%%%%%%%%%%%%%% Classification %%%%%%%%%%%%%%%%%%%%%%%%%

if visual_en;fprof = figure;end

pk_w_th = 0;
i_c0 = 1;
i_c1 = 1;

class0 =[];
class1 =[];

f = waitbar(0,'Cmod anypeak progress...');
for i_prof = 1:n_prof  % Main for loop across profiles

    pk_h_th = 10*nanstd(noise_prof(i_prof,t_basal_0:t_basal_n));
    [pks,pks_loc,pks_w,pks_prom] = findpeaks(ca_prof_smooth(i_prof,shft_l:end-shft_r),'MinPeakProminence',pk_h_th,'WidthReference','halfprom');

    if isempty(pks)
        class0(i_c0,1) = i_prof;
        i_c0 = i_c0+1;
    else
        pks = pks(pks_w>=pk_w_th);
        if isempty(pks)
            class0(i_c0,1) = i_prof;
            i_c0 = i_c0+1;
        else
            pks_loc = pks_loc(pks_w>=pk_w_th);
            pks_w = pks_w(pks_w>=pk_w_th);
        
            pks_loc = pks_loc+shft_l-1;
            class1(i_c1,1) = i_prof;
            i_c1 = i_c1+1;
        end
    end         
    if visual_en
        figure(fprof);        
        clf('reset');        
        subplot(2,1,1)        
        plot(t,ca_prof_smooth(i_prof,:),'LineWidth',1); hold on;
%         plot(t,ca_prof_base(i_prof,:),'Color',[0.75 0.75 1]);
        plot(t(pks_loc),ca_prof_smooth(i_prof,pks_loc),'*r')
        ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
        title(strcat("Prof ",num2str(i_prof)))               
        subplot(2,1,2)
        hold on;
        plot(t,ca_prof_smooth_plain(i_prof,:),'r',t,noise_prof(i_prof,:),'b');         
        plot(t(pks_loc),ca_prof_smooth_plain(i_prof,pks_loc),'*r')
        plot(t,pk_h_th*ones(length(t)),'--k')
        hold off;
        pause(0.2)          
    end
    waitbar(i_prof/n_prof,f,'Cmod anypeak progress...');
end
close(f)

if visual_en
    fclass0 = figure;
    plot(t,ca_prof_smooth(class0,:)); hold on
%     stem(t_start,max(max(ca_prof_smooth)),'r')
    title("Class no peaks")
    ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])

    fclass1=figure;
    plot(t,ca_prof_smooth(class1,:)); hold on
%     stem(t_start,max(max(ca_prof_smooth)),'r')
    title("Class with peaks")
    ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
end