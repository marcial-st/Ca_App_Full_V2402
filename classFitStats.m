function [stats] = classFitStats(profiles,filter_window,fitmodel)
% This function extracts stastistics (Mu,Sigma) from
% coefficients of fit profiles per class
% the calcium profiles are smoothed before fit
% filter_window might be compliant with the sampling period, i.e., time in
% between frames

[n_frames,n_prof] = size(profiles);
Ts = 3;
% n_params = 4;
% coeffs = zeros(n_prof,n_params);
t = (0:1:n_frames-1)*Ts;

fig_fit = figure; hold on;
xlabel('Frame'); ylabel('Ratio');
xlim([0 n_frames]);
title('Synthetic Calcium Profiles exp fit');

% f = waitbar(0,'Please wait...');

for i=1:1:n_prof    
    prof_smooth = smooth(profiles(:,i),filter_window);
    switch(fitmodel)
        case 'secord'
            fittype_secord = fittype("(1-(((exp(-z*w.*t))/sqrt(1-z^2)).*sin((w*sqrt(1-z^2)).*t+atan(sqrt(1-z^2)/z))))+b",... 
                                 dependent="y",independent="t",...
                                 coefficients=["z" "w" "b"])
            fo_secord = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0.00001,0.0001,-10],'Upper',[2,50,10]);
            cfit = fit(t',prof_smooth,fittype_secord,fo_secord);
        otherwise
            cfit = fit(t',prof_smooth,fitmodel);
    end    
    if i==1
        coeffs_aux = coeffvalues(cfit);
        n_params = length(coeffs_aux);
        coeffs = zeros(n_prof,n_params);
    end
    coeffs(i,:) = coeffvalues(cfit);
    prof_fit(i,:) = genProfileFit(coeffs(i,:),t,fitmodel,0,0);
    
    % switch(fitmodel)
    %     case 'sin3'
    %         prof_fit = coeffs(i,1)*sin(coeffs(i,2)*t+coeffs(i,3))+coeffs(i,4)*sin(coeffs(i,5)*t+coeffs(i,6))+coeffs(i,7)*sin(coeffs(i,8)*t+coeffs(i,9));
    %     case 'gauss3'
    %         prof_fit  = coeffs(i,1).*exp(-((t-coeffs(i,2))/coeffs(i,3)).^2) + coeffs(i,4).*exp(-((t-coeffs(i,5))/coeffs(i,6)).^2) + coeffs(i,7).*exp(-((t-coeffs(i,8))/coeffs(i,9)).^2);
    %     case 'poly3'
    %         prof_fit = coeffs(i,1)*t.^3+coeffs(i,2)*t.^2+coeffs(i,3)*t+coeffs(i,4);
    %     case 'poly5'
    %          prof_fit = coeffs(i,1)*t.^5+coeffs(i,2)*t.^4+coeffs(i,3)*t.^3+coeffs(i,4)*t.^2+coeffs(i,5)*t+coeffs(i,6);
    %     otherwise
    %         error('classFitStats: Fit model not found')
    % end

    prof_fit_unbias(i,:) = prof_fit(i,:)-mean(prof_fit(i,:));
    plot(t,prof_fit(i,:)+(i*0.5),'Color',[0.7 0.7 0.7]); 
    plot(t,prof_smooth+(i*0.5),'m');
    xlim([t(1) t(end)]); %ylim([0.75 1.5])
    % pause(0.05) 
    % bar_messagef = strcat('Extracting stats per class ',num2str(i));    
    % waitbar(i/n_prof,f,bar_messagef); 
end
hold off;
%% Thesis plot begin
    size_axis = 20;
    size_label = 24;
    size_legend = 14;
prof_fit_unbias_base = mean(prof_fit_unbias,1);
cfit = fit(t',prof_fit_unbias_base',fitmodel);
coeffs_base = coeffvalues(cfit);
prof_fit_unbias_base_fit = genProfileFit(coeffs_base,t,fitmodel,0,0);
figure;
hold on;
for i=1:1:n_prof
    plot(t,prof_fit_unbias,'Color',[0.75 0.75 0.75])
end
plot(t,prof_fit_unbias_base,'m','LineWidth',2)
% plot(t,prof_fit_unbias_base_fit,'c--','LineWidth',2)
ylim([-0.15 0.15])
xlim([t(1) t(end)])
             xlabel("Frame",'FontSize',size_label)
             ylabel("Magnitude",'FontSize',size_label)
set(gca,'FontSize',size_axis)
    qw{1} = plot(nan,"-",LineWidth=1,Color=[0.75 0.75 0.75]);
    qw{2} = plot(nan,"-",LineWidth=1,Color = 'm');
    legend([qw{:}],{'Unbiased fit ICP','Mean S_ICP'},'FontSize',size_legend,'Interpreter','none')  
hold off;
%% Thesis plot end
stats = coeffs_base;

% for i=1:1:n_params
%     aux = rmoutliers(coeffs(:,i));
%     stats_mean(1,i) = mean(aux);
%     stats_std(1,i) = std(aux);
% end
% stats = [stats_mean',stats_std'];
% close(f)
