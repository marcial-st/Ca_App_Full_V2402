function [class1,class0,c_vector,plateau_ratio] = cmodPlateu2(ca_prof_smooth,t_inj,t_basal_0,t_basal_n,th_plateu,visual_en,w_mainpk)
%%%%%%%%%%%%%%%%%%%%%%%%%% Vars & Params %%%%%%%%%%%%%%%%%%%%%%%%%%
% visual_en = 0;  % enable plots
[n_prof,n_frames] = size(ca_prof_smooth);
t = 1:1:n_frames;
plateau_ratio = zeros(n_prof,1);
c_vector = false(n_prof,1);
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
            muplat = nanmean(ca_prof_smooth(i_prof,t_inj+w_mainpk:end));
            max_peak = max(ca_prof_smooth(i_prof,t_inj:t_inj+w_mainpk));
            plateau_rat = ((muplat-mubase)/(max_peak-mubase));
           if plateau_rat>=th_plateu
                class1(i_c1,1) = i_prof;
                c_vector(i_prof) = true;
                if plateau_rat>1
                    plateau_ratio(i_prof) = 1;
                else
                    plateau_ratio(i_prof) = plateau_rat;
                end
                i_c1 = i_c1+1;                    
            else
                class0(i_c0,1) = i_prof;
                i_c0 = i_c0+1;                      
            end            
                                 
    waitbar(i_prof/n_prof,f,'Cmod plateu progress...');
end
close(f)

if visual_en
    fclass0=figure;
    plot(t,ca_prof_smooth(class0,:)); hold on
    xline(t_inj,'r--')
    xline(t_inj+w_mainpk,'r--')
    title("Class no plateau")
    ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])

    fclass1=figure;
    plot(t,ca_prof_smooth(class1,:)); hold on
    xline(t_inj,'r--')
    xline(t_inj+w_mainpk,'r--')
    title("Class with plateau")
    ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
end