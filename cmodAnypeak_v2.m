function [class1,class0] = cmodAnypeak_v2(profiles_matrix,visual_en,shft_l,shft_r,wfilter,ksec)
%%%%%%%%%%%%%%%%%%%%%%%%%% Vars & Params %%%%%%%%%%%%%%%%%%%%%%%%%%
% visual_en = 1;  % enable plots
fsys= readfis('fuzzy_threshold_v1.fis');
[n_prof,n_frames] = size(profiles_matrix);
t = 1:1:n_frames;
%%%%%%%%%%%%%%%%%%%%%%%%%% Profile sets %%%%%%%%%%%%%%%%%%%%%%%%%%%
ca_prof_smooth = zeros(n_prof,n_frames);                                   %% Smooth profiles with moving average filter
diff_prof_smooth = zeros(n_prof,n_frames);
ca_prof_smooth2 = zeros(n_prof,n_frames);
ca_prof_base = zeros(n_prof,n_frames);
ca_prof_smooth_plain = zeros(n_prof,n_frames);
noise_prof = zeros(n_prof,n_frames);

for i_smooth = 1:1:n_prof    
    s_prof = sgolayfilt(profiles_matrix(i_smooth,:),3,wfilter);
    noise_prof(i_smooth,:) = profiles_matrix(i_smooth,:)-s_prof;
%     i_std = nanstd(s_prof(end-20:end));
%     p_min = islocalmin(s_prof,'MinProminence',3*i_std);                                 %% find local minima 
%     p_min = islocalmin(s_prof,'MinProminence',1*nanstd(noise_prof(i_smooth,:)));        %% find local minima 
    p_min = islocalmin(s_prof);  
    p_min(1)= 1; p_min(end)= 1;                                               %% tie first/last min points 
    ca_prof_base(i_smooth,:) = interp1(t(p_min),s_prof(p_min),t,'linear'); %% Profile base, bottom envelope       
    ca_prof_smooth(i_smooth,:) = s_prof;
    diff_prof_smooth(i_smooth,2:end) = diff(s_prof);   
%     ca_prof_smooth_plain(i_smooth,:) = abs(s_prof-ca_prof_base(i_smooth,:));    %% Profile - bottom envelop to remark peaks
    ca_prof_smooth_plain(i_smooth,:) = (s_prof-ca_prof_base(i_smooth,:)).*((s_prof-ca_prof_base(i_smooth,:))>0);    %% Profile - bottom envelop to remark peaks
end

%%%%%%%%%%%%%%%%%%%%%%%% Classification %%%%%%%%%%%%%%%%%%%%%%%%%

if visual_en;fprof = figure;end

% pk_h_th = 10*mean(nanstd(ca_prof_smooth_plain(:,end-20:end),0,2));  % Peaks Threshold
% pk_w_th = 2.5;
pk_w_th = 0;
i_c0 = 1;
i_c1 = 1;
% shft = 10;
class0 =[];
class1 =[];
th1i=1;th2i=1;
f = waitbar(0,'Cmod anypeak progress...');
for i_prof = 1:n_prof  % Main for loop across profiles
%     pk_h_th = 20*nanstd(ca_prof_smooth_plain(i_prof,end-50:end));  % Peaks Threshold
%     pk_h_th = 4*nanstd(ca_prof_smooth_plain(i_prof,:));  % Peaks Threshold
%     if ((std(ca_prof_smooth(i_prof,:))/std(noise_prof(i_prof,:)))<2.5)
%         pk_h_th = 2*ksec*nanstd(noise_prof(i_prof,end-100:end));
%         th2(th2i)=i_prof;
%         th2i=th2i+1;
%     else
%         pk_h_th = ksec*nanstd(noise_prof(i_prof,end-100:end)); 
%         th1(th1i)=i_prof;
%         th1i=th1i+1;
%     end     
    [pks2,pks_loc2,pks_w2,pks_prom2] = findpeaks(ca_prof_smooth(i_prof,shft_l:end-shft_r),'WidthReference','halfprom');       
    pk_h_th = evalfis(fsys,[std(ca_prof_smooth(i_prof,shft_l:shft_l+20)) max(pks_prom2)]);
    th_vector(i_prof) = pk_h_th;
    
    [pks,pks_loc,pks_w,pks_prom] = findpeaks(ca_prof_smooth(i_prof,shft_l:end-shft_r),'MinPeakProminence',pk_h_th,'WidthReference','halfprom');    
%     [pks,pks_loc,pks_w] = findpeaks(ca_prof_smooth_plain(i_prof,shft:end),'MinPeakHeight',pk_h_th,'WidthReference','halfprom');    
%     figure;
%     plot(t,ca_prof_smooth(i_prof,:),t(pks_loc2+shft_l),ca_prof_smooth(i_prof,pks_loc2+shft_l));
%     figure;
%     plot(pks_prom2);hold on;
%     yline(pk_h_th);
%     hold off
    
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
figure; plot(th_vector,'DisplayName','Fuzzy Thresholds')
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