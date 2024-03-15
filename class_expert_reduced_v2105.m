function [pclass] = class_expert_reduced_v2105(profiles_matrix,params,visual_en)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Background %%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear all;close all;clc
% load prof_mat_3.mat
% zcIdx = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %function to return zero indices
%%%%%%%%%%%%%%%%%%%%%%%%%% Vars & Params %%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% visual_en = 1;  % enable plots
visual_sub_en = 0;  % enable plots functions

eps_val = 0.01;

% t_start = 50;   % start of experiment
% w_mainpk = 40;
% th_plateu = 0.3;
% shft = 10;

t_start = params(1);   % start of experiment
w_mainpk = params(2);
th_plateu = params(3);
e_plateu = params(4);
shft_l = params(5);
shft_r = params(6);
wfilter = params(7);
kmain = params(8);
ksec = params(9);
t_basal_0 = params(10);
t_basal_n = params(11);

cindex = [0,1,2,3,4];
clabels = ["Unknown"; %0
    "Non response"; %1
    "Basal"; %2
    "Response"; %3
    "Basal+Response";%4
    ];
cdict = containers.Map(cindex,clabels);

[n_prof,n_frames] = size(profiles_matrix);
p_idx = 1:1:n_prof;
t = 1:1:n_frames;

pclass = zeros(n_prof,1);

prof_mtx_aux = profiles_matrix;
[class0,class1,class2,class3]  = cmod_Basal_Response(prof_mtx_aux,t_start,visual_sub_en,shft_l,shft_r,wfilter,ksec,t_basal_0,t_basal_n);
pclass(p_idx(class0)) = 1;
pclass(p_idx(class1)) = 2;
pclass(p_idx(class2)) = 3;
pclass(p_idx(class3)) = 4;

% [class_basal_peak,class_no_basal_peak] = cmodBasalpeak(prof_mtx_aux(p_idx(class_with_peaks),:),t_start,visual_sub_en,shft_l,shft_r,wfilter,ksec,t_basal_0,t_basal_n);
% pclass(p_idx(class_with_peaks)) = 2; %for debug

% [class_mpeak,class_wo_mpeak]       = cmodMainpeak(prof_mtx_aux(p_idx(class_with_peaks(class_no_basal_peak)),:),t_start,w_mainpk,visual_sub_en,shft_l,shft_r,wfilter,kmain,t_basal_0,t_basal_n);
% [class_speak,class_wo_speak]       = cmodSecpeak(prof_mtx_aux(p_idx(class_with_peaks(class_no_basal_peak(class_mpeak))),:),t_start,visual_sub_en,shft_l,shft_r,wfilter,ksec,t_basal_0,t_basal_n);
% [class_op_plt,class_op_no_plt]     = cmodPlateu(prof_mtx_aux(p_idx(class_with_peaks(class_no_basal_peak(class_mpeak(class_wo_speak)))),:),t_start,t_basal_0,t_basal_n,th_plateu,e_plateu,shft_l,shft_r,visual_sub_en,w_mainpk,wfilter);
% pclass(p_idx(class_with_peaks(class_no_basal_peak(class_mpeak(class_wo_speak(class_op_no_plt)))))) = 2;
% pclass(p_idx(class_with_peaks(class_no_basal_peak(class_mpeak(class_wo_speak(class_op_plt)))))) = 3;
% 
% [class_op_plt,class_op_no_plt]     = cmodPlateu(prof_mtx_aux(p_idx(class_with_peaks(class_no_basal_peak(class_mpeak(class_speak)))),:),t_start,t_basal_0,t_basal_n,th_plateu,e_plateu,shft_l,shft_r,visual_sub_en,w_mainpk,wfilter);
% pclass(p_idx(class_with_peaks(class_no_basal_peak(class_mpeak(class_speak(class_op_no_plt)))))) = 4;
% pclass(p_idx(class_with_peaks(class_no_basal_peak(class_mpeak(class_speak(class_op_plt)))))) = 5;
% 
% [class_speak,class_wo_speak]       = cmodSecpeak(prof_mtx_aux(p_idx(class_with_peaks(class_no_basal_peak(class_wo_mpeak))),:),t_start,visual_sub_en,shft_l,shft_r,wfilter,ksec,t_basal_0,t_basal_n);
% [class_op_plt,class_op_no_plt]     = cmodPlateu(prof_mtx_aux(p_idx(class_with_peaks(class_no_basal_peak(class_wo_mpeak(class_speak)))),:),t_start,t_basal_0,t_basal_n,th_plateu,e_plateu,shft_l,shft_r,visual_sub_en,w_mainpk,wfilter);
% pclass(p_idx(class_with_peaks(class_no_basal_peak(class_wo_mpeak(class_speak(class_op_no_plt)))))) = 6;
% pclass(p_idx(class_with_peaks(class_no_basal_peak(class_wo_mpeak(class_speak(class_op_plt)))))) = 7;
% 
% [class_mpeak,class_wo_mpeak]       = cmodMainpeak(prof_mtx_aux(p_idx(class_with_peaks(class_basal_peak)),:),t_start,w_mainpk,visual_sub_en,shft_l,shft_r,wfilter,kmain,t_basal_0,t_basal_n);
% [class_speak,class_wo_speak]       = cmodSecpeak(prof_mtx_aux(p_idx(class_with_peaks(class_basal_peak(class_mpeak))),:),t_start,visual_sub_en,shft_l,shft_r,wfilter,ksec,t_basal_0,t_basal_n);
% [class_op_plt,class_op_no_plt]     = cmodPlateu(prof_mtx_aux(p_idx(class_with_peaks(class_basal_peak(class_mpeak(class_wo_speak)))),:),t_start,t_basal_0,t_basal_n,th_plateu,e_plateu,shft_l,shft_r,visual_sub_en,w_mainpk,wfilter);
% pclass(p_idx(class_with_peaks(class_basal_peak(class_mpeak(class_wo_speak(class_op_no_plt)))))) = 8;
% pclass(p_idx(class_with_peaks(class_basal_peak(class_mpeak(class_wo_speak(class_op_plt)))))) = 9;
% 
% [class_op_plt,class_op_no_plt]     = cmodPlateu(prof_mtx_aux(p_idx(class_with_peaks(class_basal_peak(class_mpeak(class_speak)))),:),t_start,t_basal_0,t_basal_n,th_plateu,e_plateu,shft_l,shft_r,visual_sub_en,w_mainpk,wfilter);
% pclass(p_idx(class_with_peaks(class_basal_peak(class_mpeak(class_speak(class_op_no_plt)))))) = 10;
% pclass(p_idx(class_with_peaks(class_basal_peak(class_mpeak(class_speak(class_op_plt)))))) = 11;
% 
% [class_speak,class_wo_speak]       = cmodSecpeak(prof_mtx_aux(p_idx(class_with_peaks(class_basal_peak(class_wo_mpeak))),:),t_start,visual_sub_en,shft_l,shft_r,wfilter,ksec,t_basal_0,t_basal_n);
% [class_op_plt,class_op_no_plt]     = cmodPlateu(prof_mtx_aux(p_idx(class_with_peaks(class_basal_peak(class_wo_mpeak(class_speak)))),:),t_start,t_basal_0,t_basal_n,th_plateu,e_plateu,shft_l,shft_r,visual_sub_en,w_mainpk,wfilter);
% pclass(p_idx(class_with_peaks(class_basal_peak(class_wo_mpeak(class_speak(class_op_no_plt)))))) = 12;
% pclass(p_idx(class_with_peaks(class_basal_peak(class_wo_mpeak(class_speak(class_op_plt)))))) = 13;
% 
% [class_op_plt,class_op_no_plt]     = cmodPlateu(prof_mtx_aux(p_idx(class_with_peaks(class_basal_peak(class_wo_mpeak(class_wo_speak)))),:),t_start,t_basal_0,t_basal_n,th_plateu,e_plateu,shft_l,shft_r,visual_sub_en,w_mainpk,wfilter);
% pclass(p_idx(class_with_peaks(class_basal_peak(class_wo_mpeak(class_wo_speak(class_op_no_plt)))))) = 14;
% pclass(p_idx(class_with_peaks(class_basal_peak(class_wo_mpeak(class_wo_speak(class_op_plt)))))) = 15;
% 
% pclass(pclass==0)= 6; % Remaining ones have only sec peaks 

pclass_u = unique(pclass);
if visual_en
    for i=1:length(pclass_u)
        figure;
        plot(t,profiles_matrix(p_idx(pclass == pclass_u(i)),:))
        title(strcat("Class ",num2str(pclass_u(i)),": ",values(cdict,{pclass_u(i)})))
    end
end

for i=1:1:length(pclass_u)
    NCLS(i).cVector = find(pclass==pclass_u(i))';
    string_class(i) = {strcat("Class ",num2str(pclass_u(i)))};    
end 
disp('Class Expert Finished!')