function [ca_prof_smooth,noise_prof] = smoothProfiles(profiles_matrix,wfilter)

[n_prof,n_frames] = size(profiles_matrix);
ca_prof_smooth = zeros(n_prof,n_frames); 
noise_prof = zeros(n_prof,n_frames); 

for i_smooth = 1:1:n_prof    
    s_prof = sgolayfilt(profiles_matrix(i_smooth,:),3,wfilter);
    noise_prof(i_smooth,:) = profiles_matrix(i_smooth,:)-s_prof;
    ca_prof_smooth(i_smooth,:) = s_prof;
end