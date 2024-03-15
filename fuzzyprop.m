% clear;clc;
load data_prop.mat

ca_prof_smooth(:,1:10) = [];
diff_prof_smooth(:,1:10) = [];
noise_prof(:,1:10) = [];
t(end-9:end) = [];

[nprof,nframe]=size(ca_prof_smooth);

figure; plot(t,ca_prof_smooth); title('Ca prof smooth'); xlim([1 length(t)])
figure; plot(t,diff_prof_smooth); title('Diff prof smooth'); xlim([1 length(t)])
figure; plot(t,noise_prof); title('Noise prof'); xlim([1 length(t)])


% f1=figure;
% for i=1:1:nprof
%     hold on;
%     plot(t,ca_prof_smooth(i,:)-mean(ca_prof_smooth(i,:)),'DisplayName','Ca prof'); 
%     plot(t,diff_prof_smooth(i,:),'DisplayName','Diff prof'); xlim([1 length(t)])
%     plot(t,noise_prof(i,:),'DisplayName','Noise prof'); xlim([1 length(t)])
%     legend
%     xlim([1 length(t)]);ylim([-0.2 0.2]);
%     pause()
%     clf(f1)
% end

mu_prof = nanmean(ca_prof_smooth,2);
mu_diff = nanmean(diff_prof_smooth,2);
mu_noise = nanmean(noise_prof,2);

figure; 
subplot(3,1,1)
plot(mu_prof);
title('Mean prof')
subplot(3,1,2)
plot(mu_diff);
title('Mean diff')
subplot(3,1,3)
plot(mu_noise);
title('Mean noise')

sigma_prof = std(ca_prof_smooth,0,2);
sigma_diff = std(diff_prof_smooth,0,2);
sigma_noise = std(noise_prof,0,2);

figure; 
subplot(3,1,1)
plot(sigma_prof);
title('Std prof')
subplot(3,1,2)
plot(sigma_diff);
title('Std diff')
subplot(3,1,3)
plot(sigma_noise);
title('Std noise')

max_prof = max(ca_prof_smooth,[],2);
max_diff = max(diff_prof_smooth,[],2);
max_noise = max(noise_prof,[],2);

figure; 
subplot(3,1,1)
plot(max_prof);
title('Max prof')
subplot(3,1,2)
plot(max_diff);
title('Max diff')
subplot(3,1,3)
plot(max_noise);
title('Max noise')


min_prof = min(ca_prof_smooth,[],2);
min_diff = min(diff_prof_smooth,[],2);
min_noise = min(noise_prof,[],2);

figure; 
subplot(3,1,1)
plot(min_prof);
title('Min prof')
subplot(3,1,2)
plot(min_diff);
title('Min diff')
subplot(3,1,3)
plot(min_noise);
title('Min noise')

snr = sigma_prof./sigma_noise;
figure;
plot(snr)
title('SNR')

v_prom = zeros(1,nprof);
for i_prof =1:1:nprof
   [pks,pks_loc,pks_w,pks_prom] = findpeaks(ca_prof_smooth(i_prof,:),'WidthReference','halfprom');
   v_prom(i_prof)= max(pks_prom);
end

figure;
plot(v_prom)
title('Prominence')

% fsys= readfis('fuzzy_threshold_v1.fis');
% pk_h_th = evalfis(fsys,[std(diff_prof_smooth(i_prof,:)) max(diff_prof_smooth(i_prof,:))]);