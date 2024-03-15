% function [CELLSf] = featureExtraction_v200927(profiles_matrix)
clear all;close all;clc
load prof_mat_1.mat
% pr_mx=profiles_matrix';
visual_en = 0;
[n_prof,n_frames] = size(profiles_matrix);
t = 1:1:n_frames;
idx_vector = 1:1:n_prof;



for i_smooth = 1:1:n_prof
    s_prof = smooth(profiles_matrix(i_smooth,:),5)';      
    ca_prof_smooth(i_smooth,:) = s_prof;    
end

pr_mx=ca_prof_smooth;


corr_th = 0.95;



i_c = 1;
while ((length(idx_vector)>1)&&(i_c<30))
    corr_mat = corrcoef(pr_mx);
    good_ones = corr_mat(:,1)>=corr_th;
    class_i=idx_vector(good_ones);
    
    NCLS(i_c).vector=class_i;
    
    idx_vector(good_ones)= [];    
    pr_mx(:,good_ones)=[];
%     good_ones = [];
%     corr_mat = [];
    
    figure;plot(t,profiles_matrix(class_i,:));title(strcat("Class ",num2str(i_c)));
    
    i_c = i_c+1;
end
