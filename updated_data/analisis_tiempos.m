clear all;close all;clc

path_main="C:\Users\marci\OneDrive\Documentos\MATLAB\sandbox\Ca_App_V2-main\updated_data";
exp_name=["100528_1.mat";
          "100528_4.mat";
          "100528_5.mat";
          "100528_6.mat";
          "100528_8.mat";
          "100601_1.mat"];

Ts = 3;
delta_frames = 0;

cells_global = [];
v_fstim = zeros(length(exp_name),1);
v_fstart = zeros(length(exp_name),1);
v_fend = zeros(length(exp_name),1);
v_nprof = zeros(length(exp_name),1);

for i_exp = 1:1:length(exp_name)
    path_exp = fullfile(path_main,exp_name(i_exp));
    load(path_exp);        

    v_fstim(i_exp) = floor(exp_data.params.tstim/Ts);
    v_fstart(i_exp) = floor(exp_data.params.tstart/Ts);
    v_fend(i_exp) = floor(exp_data.params.tend/Ts);    
end

fstim_ref = min(v_fstim);
fend_ref = min(v_fend-v_fstim);
t = [0:1:(fstim_ref+fend_ref)-1]*Ts;

figure; hold on;
ca_prof_global = [];
for i_exp = 1:1:length(exp_name)
% for i_exp = 1
    path_exp = fullfile(path_main,exp_name(i_exp));
    load(path_exp);    
    [cells_available,~] = get_offline_cells(delta_frames,exp_data,Ts);       
    ca_prof_i = reshape([cells_available.ca_profile]',length(cells_available(1).ca_profile),length(cells_available))';    
    ca_prof_i = ca_prof_i(:,v_fstim(i_exp)-fstim_ref+1:v_fstim(i_exp)+fend_ref);
    ca_prof_global = [ca_prof_global;ca_prof_i];
    plot(t,ca_prof_i+(1*(i_exp-1)));
end
xline(fstim_ref*Ts,'k--')
xlim([t(1) t(end)]);
title("All Ca profiles aligned+offset")
hold off;

figure;
plot(t,ca_prof_global)
title("All Ca profiles aligned")
xlim([t(1) t(end)]);