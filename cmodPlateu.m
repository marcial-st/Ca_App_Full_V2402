function [class1,class0] = cmodPlateu(profiles_matrix,t_start,t_basal_0,t_basal_n,th_plateu,e_plateu,shft_l,shft_r,visual_en,w_mainpk,wfilter)
%%%%%%%%%%%%%%%%%%%%%%%%%% Vars & Params %%%%%%%%%%%%%%%%%%%%%%%%%%
% visual_en = 0;  % enable plots
[n_prof,n_frames] = size(profiles_matrix);
t = 1:1:n_frames;
%%%%%%%%%%%%%%%%%%%%%%%%%% Profile sets %%%%%%%%%%%%%%%%%%%%%%%%%%%

[ca_prof_smooth,noise_prof] = smoothProfiles(profiles_matrix,wfilter);

%%%%%%%%%%%%%%%%%%%%%%%% Classification %%%%%%%%%%%%%%%%%%%%%%%%%

% if visual_en;fprof = figure;end

% pk_h_th = 10*mean(nanstd(ca_prof_smooth_plain(:,end-20:end),0,2));  % Peaks Threshold
% pk_h_th = 0.035
pk_w_th = 0;
i_c0 = 1;
i_c1 = 1;
% shft = 10;
class0 = [];
class1 = [];
f = waitbar(0,'Cmod plateu progress...');
for i_prof = 1:n_prof  % Main for loop across profiles
         
            mubase = nanmean(ca_prof_smooth(i_prof,t_basal_0:t_basal_n));
            muplat = nanmean(ca_prof_smooth(i_prof,t_start+w_mainpk:end-shft_r));
            
           if (muplat-mubase)>=th_plateu
                class1(i_c1,1) = i_prof;
                i_c1 = i_c1+1;                    
            else
                class0(i_c0,1) = i_prof;
                i_c0 = i_c0+1;                      
            end            
            
            
            
%             p0 = mean(ca_prof_smooth(i_prof,t_start-10:t_start));
% %             p1 = mean(ca_prof_smooth(i_prof,t_start+w_mainpk:t_start+w_mainpk+20));
% %             p2 = mean(ca_prof_smooth(i_prof,end-50:end));
%             p1 = mean(ca_prof_smooth(i_prof,t_start+w_mainpk:t_start+w_mainpk+10));
%             p2 = mean(ca_prof_smooth(i_prof,end-10:end));
%             
%             if (((p1-p0)>th_plateu) && ((p2-p0)>th_plateu))
%                 if abs(p1-p2)<=e_plateu
%                     class1(i_c1,1) = i_prof;
%                     i_c1 = i_c1+1;                    
%                 else
%                     class0(i_c0,1) = i_prof;
%                     i_c0 = i_c0+1;                      
%                 end
%             else
%                 class0(i_c0,1) = i_prof;
%                 i_c0 = i_c0+1;  
%             end
                       
            
            
% %             if abs(p2-p1)>=(th_plateu*(ca_prof_smooth(i_prof,(pks_loc(1)))-p1))
%             if p2-p0>=th_plateu
%                 class1(i_c1,1) = i_prof;
%                 i_c1 = i_c1+1;
%             else
%                 class0(i_c0,1) = i_prof;
%                 i_c0 = i_c0+1;            
%             end
% % %         else
% % %             class0(i_c0,1) = i_prof;
% % %             i_c0 = i_c0+1;            
% % %         end 
% % %     else
% % %         class0(i_c0,1) = i_prof;
% % %         i_c0 = i_c0+1;
% % %     end    
% %     if visual_en
% %         figure(fprof);        
% %         clf('reset');        
% %         subplot(1,2,1)
% %         hold on;
% %         plot(t,ca_prof_smooth_plain(i_prof,:));         
% %         plot(t(pks_loc),ca_prof_smooth_plain(i_prof,pks_loc),'*r')
% %         hold off;
% %         subplot(1,2,2)
% %         plot(t,ca_prof_smooth(i_prof,:),'LineWidth',1); hold on;
% %         plot(t,ca_prof_base(i_prof,:),'Color',[0.75 0.75 1]);
% %         plot(t(pks_loc),ca_prof_smooth(i_prof,pks_loc),'*r')
% %         ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
% %         title(strcat("Prof ",num2str(i_prof)))
% %         pause(0.2)          
% %     end
    waitbar(i_prof/n_prof,f,'Cmod plateu progress...');
end
close(f)

if visual_en
    fclass0=figure;
    plot(t,ca_prof_smooth(class0,:)); hold on
    stem(t_start,max(max(ca_prof_smooth)),'r')
    title("Class no plateau")
    ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])

    fclass1=figure;
    plot(t,ca_prof_smooth(class1,:)); hold on
    stem(t_start,max(max(ca_prof_smooth)),'r')
    title("Class with plateau")
    ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
end
