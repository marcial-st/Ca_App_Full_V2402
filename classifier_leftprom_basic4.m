function [pclass,global_events_bin,global_events] = classifier_leftprom_basic4(profiles_smooth,th_artifact,n_frames,Ts,th_diff,params,th_prom,peak_width_f,en_artifact_detector_2)
en_art_noreg2 = 1;
en_shift_8020 = 1; th_shift_p=0.20;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-----------      Loading Profiles      -----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Start: Loading profiles")            


[n_prof,n_frames] = size(profiles_smooth);
peak_width_vector = zeros(1,n_prof);
t = Ts.*linspace(0,n_frames-1,n_frames);
class_vector = string(zeros(n_prof,2));                    
  
profiles_diff = diff(profiles_smooth,1,2)>th_diff;

if isfield(params,'tstim') && ~isempty(params.tstim)
    inj_time  = params.tstim;
    inj_frame = floor(inj_time/3);
    disp("     Injury time got from params")
else
    [~,inj_frame] = max(sum(profiles_diff,1));
    inj_time = t(inj_frame)-36; %45,36

    inj_frame = floor(inj_time/3);
    disp("     Injury time calculated")
end

global_events_val = zeros(n_prof,n_frames);
global_events = zeros(n_prof,n_frames);
global_events_w = zeros(n_prof,n_frames);
global_events_tri = zeros(n_prof,n_frames);

global_events_inv = zeros(n_prof,n_frames);
global_events_w_inv = zeros(n_prof,n_frames);
global_events_tri_inv = zeros(n_prof,n_frames);

basal_line = zeros(n_prof,1);

disp("End: Loading profiles")				
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%----------- Photobleaching compensation -----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
disp("Start: Photobleaching compensation")
[std_prof, std_idx]=sort(std(profiles_smooth,0,2),'ascend');
n_comp = 20;
l_fit = zeros(n_comp,2);
for i_comp = 1:n_comp            
    l_fit(i_comp,:) = polyfit(t,profiles_smooth(std_idx(i_comp),:),1);
    yfit = l_fit(i_comp,1)*t+l_fit(i_comp,2);
end
l_fit = mean(l_fit);
y_comp = (l_fit(1)*t);
profiles_smooth_comp = profiles_smooth-(ones(n_prof,1)*y_comp);
profiles_smooth_bu = profiles_smooth;
profiles_smooth = profiles_smooth_comp;
disp("End: Photobleaching compensation")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%----------- Prominence thresholding -----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_stop = floor((params.tend)/Ts);
if f_stop>n_frames
    f_stop = n_frames;
    disp("Warning! f_stop > n_frames, then f_stop = n_frames was assigned")
end
disp("Start: Prominence thresholding")
for i_prof=1:n_prof
    [pks_all,pk_max_loc_f,pks_w_all,pks_prom_all,tri_vector] = findpeaks_wrapper(profiles_smooth(i_prof,1:f_stop),0,i_prof);        
%     [pks_all_inv,pk_max_loc_f_inv,pks_w_all_inv,pks_prom_all_inv,tri_vector_inv] = findpeaks_wrapper_inv(profiles_smooth(i_prof,1:f_stop),0,i_prof);

    global_events(i_prof,pk_max_loc_f) = pks_prom_all;
    global_events_w(i_prof,pk_max_loc_f) = pks_w_all;
    global_events_tri(i_prof,1:f_stop) = tri_vector;
    global_events_val(i_prof,pk_max_loc_f) = pks_all;

%     global_events_inv(i_prof,pk_max_loc_f_inv) = pks_prom_all_inv;
%     global_events_w_inv(i_prof,pk_max_loc_f_inv) = pks_w_all_inv;
%     global_events_tri_inv(i_prof,1:f_stop) = tri_vector_inv; 
end             

basal_line = mean(profiles_smooth(:,1:inj_frame),2);  % calculating basal levels  

% Thresholding: calcium true event
global_events_bin = global_events>th_prom;
global_events_bin_inv = global_events_inv>th_prom;

disp("Done: Prominence thresholding")
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-----------    Artifact Detection    --------NEW based on basal levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
art_approac2_en = en_artifact_detector_2
if art_approac2_en == 1
    for i_prof=1:n_prof
        mask = global_events_val(i_prof,inj_frame:inj_frame+peak_width_f)>basal_line(i_prof);
        
        global_events(i_prof,inj_frame:inj_frame+peak_width_f) = global_events(i_prof,inj_frame:inj_frame+peak_width_f).*mask;
        global_events_w(i_prof,inj_frame:inj_frame+peak_width_f) = global_events_w(i_prof,inj_frame:inj_frame+peak_width_f).*mask;
        global_events_val(i_prof,inj_frame:inj_frame+peak_width_f) = global_events_val(i_prof,inj_frame:inj_frame+peak_width_f).*mask;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-----------    Artifact Detection    -------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
disp("Start: Artifact detection")         
global_events_vector = sum(global_events_bin,1); 
global_events_vector(inj_frame:inj_frame+peak_width_f) = 0;
artifacts_vector = find(global_events_vector > th_artifact);
artifact_vector_alg =[];
if isempty(artifacts_vector)
    disp(strcat("     No artifacts found with threshold = ",num2str(th_artifact)))
else
    artifact_bin = logical(zeros(size(global_events_bin)));
    dav = diff(artifacts_vector);
    art_g_idx = find(dav>5);
    art_g_idx = [0,art_g_idx,length(artifacts_vector)];
    
    for i_group = 1:1:length(art_g_idx)-1
        art_group(1) = artifacts_vector(art_g_idx(i_group)+1);
        art_group(2) = artifacts_vector(art_g_idx(i_group+1));                        
        for i_prof = 1:1:n_prof
            mins = find(islocalmin(profiles_smooth(i_prof,:)));
            aux = sort(mins(mins<art_group(1)),'descend');
            if isempty(aux)
                art_l = art_group(1);
            elseif (en_art_noreg2 && (inj_frame<aux(1)) && (aux(1)<(inj_frame+peak_width_f)))
                art_l = inj_frame+peak_width_f;
            else
                art_l = aux(1);
            end
            
            aux = sort(mins(mins>art_group(2)),'ascend');                                                
            if isempty(aux)
                art_r = art_group(2);
            elseif (en_art_noreg2 && (inj_frame<aux(1)) && (aux(1)<(inj_frame+peak_width_f)))
                art_r = inj_frame;
            else
                art_r = aux(1);
            end
            
              art_l = art_l - 15;
              art_l = art_r + 15;

            if ((max(global_events(i_prof,art_l:art_r))<3*th_prom)) %% OR (1.5*th_prom menor que la basal)
                artifact_bin(i_prof,art_l:art_r) = logical(global_events_bin(i_prof,art_l:art_r));
                global_events_bin(i_prof,art_l:art_r) = 0;
                if not(isempty(find(artifact_bin(i_prof,art_l:art_r), 1)))
                    artifact_vector_alg =[artifact_vector_alg i_prof];
                end
            end 
        end	
    end
end
disp("End: Artifact detection")    		   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-----------    Classification        -------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp("Start: Classification")
% class_vector = string(zeros(n_prof,2));
class_vector = zeros(n_prof,2);

% spks_bin = zeros(length(exp_data.spks),n_frames);
% for k=1:length(exp_data.spks)                
    % if ~isempty(exp_data.spks(k).pks)
        % for k2=1:length(exp_data.spks(k).pks(:,1))
            % if ~(exp_data.spks(k).pks(k2,1)<=0)
                % aux_frame = floor(exp_data.spks(k).pks(k2,1)/Ts); 
                % spks_bin(k,aux_frame) = exp_data.spks(k).pks(k2,2);
            % end
        % end
    % end
% end
% spks_bin = spks_bin>0;
  % prof_resp = global_events_bin(find(gt_data(:,2)=="R"),:);
% prof_resp = spks_bin(gt_data(:,2)=="R",:);
% prof_resp(:,1:inj_frame)=0;
% [n_profresp,~] = size(prof_resp);
% first_resp = zeros(n_profresp,1);
% prof_resp=[]; % debug wa gt data for peak_width calculation
% if ~isempty(prof_resp)
    % for i_prof =1:1:n_profresp
        % aux = find(prof_resp(i_prof,:));
          % first_resp(i_prof) = aux(1);
        % if ~isempty(aux)
            % first_resp(i_prof) = aux(1);
        % else
            % first_resp(i_prof) = first_resp(i_prof-1);
            % disp("Warning: first response in xRx not found, it must be there :S")
        % end                    
    % end
    % [first_resp,~] = rmoutliers(first_resp,'mean');
    % peak_width_f = max(first_resp)+3-inj_frame;
    % peak_width = peak_width_f*Ts;   
% else
    % disp('Responding profiles not found, then fixed region 2')
    % peak_width_f = peak_width_f_bu;
    % peak_width = peak_width_f*Ts; 
% end

    peak_width = peak_width_f*Ts; 

for i_prof=1:n_prof
    
    %% en_shift_8020
    if i_prof == 1
        global_events_bin_bu = global_events_bin;
        global_events_bu     = global_events;
        global_events = global_events.*global_events_bin_bu;
    end
    idx_events = find(global_events_w(i_prof,inj_frame:inj_frame+peak_width_f));
    if ~isempty(idx_events)
        left_frame = global_events_w(i_prof,idx_events(1)+inj_frame-1)/Ts;                        
        [aux_bin_update,aux_prom_update]= UpdateGlobalEvents_th(profiles_smooth(i_prof,inj_frame-left_frame:inj_frame+peak_width_f),[zeros(1,left_frame),global_events_bin(i_prof,inj_frame:inj_frame+peak_width_f)],[zeros(1,left_frame),global_events(i_prof,inj_frame:inj_frame+peak_width_f)],[zeros(1,left_frame),global_events_w(i_prof,inj_frame:inj_frame+peak_width_f)],Ts,th_shift_p);
     if isempty(find(aux_bin_update(1:left_frame)))
         global_events_bin(i_prof,inj_frame:inj_frame+peak_width_f)=aux_bin_update(left_frame+1:end);
         global_events(i_prof,inj_frame:inj_frame+peak_width_f)=aux_prom_update(left_frame+1:end);
     else
         global_events_bin(i_prof,inj_frame-left_frame:inj_frame+peak_width_f)=aux_bin_update;
         global_events(i_prof,inj_frame-left_frame:inj_frame+peak_width_f)=aux_prom_update;
     end
    else
        disp("Warning: No activity detected in this profile")
    end
    %% end
      
    % case("basic4")
    if sum(global_events_bin(i_prof,1:inj_frame))>=1
%         class_vector(i_prof,1) = "R";
          class_vector(i_prof,1) = 1;
    else
%         class_vector(i_prof,1) = "NR";
        class_vector(i_prof,1) = 0;
    end

    if sum(global_events_bin(i_prof,inj_frame+1:n_frames))>=1
%         class_vector(i_prof,2) = "R";
        class_vector(i_prof,2) = 1;
    else
%         class_vector(i_prof,2) = "NR";
        class_vector(i_prof,2) = 0;
    end  
 
end
pclass = (class_vector(:,1)*2^0)+(class_vector(:,2)*2^1);
disp("End: Classification")