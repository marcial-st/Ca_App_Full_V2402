path = 'C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Peak_Analysis_Ground_Truth\Smooth';
files = dir(path);

nf = length(files);
wfilter = 11;

for i=1:1:nf
    file_name = files(i).name;
    if contains(file_name,'.mat','IgnoreCase',true) 
        disp(strcat("Processing file ",file_name,"   ",num2str(round(i/nf*100),"%")))
        load(fullfile(path,file_name))
        [smooth_profiles_matrix,noise_profiles_matrix] = smoothProfiles(profiles_matrix,wfilter);        
        save(fullfile(path,file_name),'profiles_matrix','smooth_profiles_matrix','noise_profiles_matrix','spks')
    end   
end