function [TP,TN,FP,FN,total,fullrates] = peakassessment(spks,profiles_smooth,profiles_original,th_p,th_w,th_k,epsilon,params,approach_sel)

% epsilon = 3;
Ts = 3;
[np,nf] = size(profiles_smooth);

%% get time params
if isfield(params,'tstart')
    if (params.tstart==0)
        tstart = 1;
    else
        tstart = params.tstart;
    end
else
    tstart = 1;
end

if isfield(params,'tstim')
    tstim = params.tstim;
else
    tstim = 20;
end

if isfield(params,'tstop')
    tstop = params.tstop;
else
    tstop = nf*Ts;
end

%% main analysis loop
gtr_i = 1;
TP = 0;
TN = 0;
FP = 0;
FN = 0;                        
total = 0;
total_gt = 0;
fullrates.TP = [];
fullrates.TN = [];
fullrates.FP = [];
fullrates.FN = [];
fullrates.total = [];
fullrates.total_g = [];
fullrates.count_pos = [];
fullrates.count_neg = [];
fullrates.filt_a_p = [];
fullrates.filt_a_n = [];
% approach_sel = 1;

tstim
tstop
tstop_f = floor(tstop/Ts)
tstim_f = floor(tstim/Ts);

% a=[abs(profiles_smooth(:,tstart:tstim_f)-profiles_original(:,tstart:tstim_f))];
% stdth = min(nanstd(a,0,2));
aux_matrix = zeros(np,tstop_f-tstart+1);

for k = 1:1:np
% for k = 90:1:90
    % two sets of peaks, auto & ground truth
    pks_gt = spks(k).pks;
    switch approach_sel
        case 2
            [pks_all,pks_loc_all,aux_all_p,aux_all_n,filt_a_p,filt_a_n]=th_fine_approach_v1(profiles_smooth(k,tstart:tstop_f),th_p);
            fullrates.filt_a_p{k} = filt_a_p;
            fullrates.filt_a_n{k} = filt_a_n;
            fullrates.pks_all{k} = pks_all;
            fullrates.pks_loc_all{k} = pks_loc_all;
            aux_all_p = aux_all_p*Ts;
            aux_all_n = aux_all_n*Ts;
        otherwise
            [pks_all,pk_max_loc_f,pks_w_all,pks_prom_all] = findpeaks_wrapper(profiles_smooth(k,tstart:tstop_f));
                % Thresholding cuantification, get logical vectors
            if not(isempty(pk_max_loc_f))
                filt_wp = (pks_w_all*Ts)>=th_w;
                filt_pp = pks_prom_all>=th_p;
%                 filt_kp = pks_prom_all>=(th_k*nanstd(profiles_smooth(k,tstart:tstim)));
                filt_kp = pks_prom_all>=(th_k*nanstd(abs(profiles_smooth(k,tstart:tstim_f)-profiles_original(k,tstart:tstim_f))));
%                 filt_kp = pks_prom_all>=(th_k*stdth);
                

                filt_a_p = logical(filt_wp.*filt_pp.*filt_kp);
                filt_a_n = not(filt_a_p);

                aux_all_p = (pk_max_loc_f(filt_a_p)*Ts)-Ts;
                aux_all_n = (pk_max_loc_f(filt_a_n)*Ts)-Ts;

                fullrates.filt_a_p{k} = filt_a_p;
                fullrates.filt_a_n{k} = filt_a_n;
                fullrates.pks_all{k} = pks_all;
                fullrates.pks_loc_all{k} = pk_max_loc_f;   
                if sum(filt_a_p)>0
                    aux_matrix(k,pk_max_loc_f(filt_a_p)) = pks_all(filt_a_p);
                end
            end   
    end 
 
    % Ground Truth cuantification, get logical vectors
    if not(isempty(pks_gt))
        n_gt = length(pks_gt(:,1));
        aux_gt = nonzeros(pks_gt(:,1));
    end
    
    % there are four cases based on the peaks availability in each set
                 % GroundTruth           %AutoDetection
    aux_e = strcat(num2str(not(isempty(pks_gt))),num2str(not(isempty(pks_all))));    
    % TP: True Positive
    % TN: True Negative
    % FP: False Positive
    % FN: False Negative
    switch aux_e
        case '00'
%             disp(strcat("No peaks neither in GroudTruth nor in AutoDetection, profile:",num2str(k)))
            TP_i = 0;
            TN_i = 0;
            FP_i = 0;
            FN_i = 0;                        
            total_i = 0;
            total_g = 0;
            count_p = 0;
            count_n = 0;
        case '01'
%             disp(strcat("No peaks in GroundTruth set, profile: ",num2str(k)))            
            TP_i = 0;
            FN_i = 0;
            TN_i = sum(filt_a_n);
            FP_i = sum(filt_a_p);
            total_i = length(pks_all);
            total_g = 0;
            count_p = sum(filt_a_p);
            count_n = sum(filt_a_n);
        case '10'
%             disp(strcat("No peaks in Autodetection set, profile: ",num2str(k)))            
            TP_i = 0;
            TN_i = 0;            
            FP_i = 0;             
            FN_i = n_gt;
            total_i = 0;
            total_g = n_gt;
            count_p = 0;
            count_n = 0;            
        case '11'
%             disp('Caso 1 1')
            match_count_p = 0;
            match_count_n = 0;
            n_match_count_n = 0;
            for j=1:1:length(aux_gt) % Los puntos de GT que no pasan el umbral se envían directamente a FN
                
                aux_match_p = (sum(abs(aux_all_p-aux_gt(j))<=epsilon));
                aux_match_n = (sum(abs(aux_all_n-aux_gt(j))<=epsilon));
                
                if (aux_match_p>0)
                    match_count_p = match_count_p+1;
                end                
                if ((aux_match_p==0)&&(aux_match_n>0))
                    match_count_n = match_count_n+1;
                end
                if ((aux_match_p==0)&&(aux_match_n==0))
                    n_match_count_n = n_match_count_n+1;
                end
            end
            total_i = length(pks_all);
            total_g = n_gt;
            TP_i = match_count_p;
            FP_i = length(aux_all_p)-match_count_p;
            if FP_i<0
                disp('FP_i<0')
            end
            TN_i = length(aux_all_n)-match_count_n;
%             FN_i = match_count_n+n_match_count_n;            
            FN_i = match_count_n;    
            if FN_i<0
                disp('FN_i<0')
            end
            count_p = sum(filt_a_p);
            count_n = sum(filt_a_n);
            GT_match_rate_i = 100*(total_g-n_match_count_n)/total_g;
            fullrates.gt_match_rate(gtr_i) = GT_match_rate_i;
            gtr_i = gtr_i+1;
            if (TP_i+TN_i+FP_i+FN_i~=total_i)
                disp(strcat('Warning: stats did not match with total_i, profile:',num2str(k)))
            end
        otherwise
            disp('This case should not be :O')
    end
    
    fullrates.TP(k) = TP_i;
    fullrates.TN(k) = TN_i;
    fullrates.FP(k) = FP_i;
    fullrates.FN(k) = FN_i;
    fullrates.total(k) = total_i;
    fullrates.total_g(k) = total_g;
    fullrates.count_neg(k) = count_n;
    fullrates.count_pos(k) = count_p;
    fullrates.aux_matrix = aux_matrix;
    
    TP = TP+TP_i;
    TN = TN+TN_i;
    FP = FP+FP_i;
    FN = FN+FN_i;                        
    total = total+total_i;
    total_gt = total_gt+total_g;
end
disp("Peak assessment done!")