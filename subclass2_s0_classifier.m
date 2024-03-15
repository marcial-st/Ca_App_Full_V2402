function [NCSL2_0] = subclass2_s0_classifier(spks,idx_vector,Ts,tstim,th_vals,profiles,exp_name,global_events,global_events_bin,en_ignore_low_resp)

colors = [0.6751    0.0851    0.2852;
          0.7538    0.5855    0.9400;
          0.1859    0.7648    0.9213;
          0.5342    0.2018    0.9242;
          0.3024    0.5672    0.1932;
          0.7793    0.8332    0.8613;
          0.2339    0.9671    0.5985;
          0.9241    0.4592    0.6695;
          0.8074    0.9329    0.1014;
          0.4725    0.2987    0.9894;
          0.8829    0.8844    0.2303;
          0.6845    0.5768    0.3458;
          0.5100    0.5680    0.2144;
          0.1551    0.3441    0.3465;
          0.6625    0.9827    0.8964;
          0.5247    0.8909    0.4667];

en_art_3 = en_ignore_low_resp;
cindex = [0,1];
clabels = ["RLow"; %0
            "RHigh"];%1
    
cdict = containers.Map(cindex,clabels);

[n_prof,n_frame]=size(profiles);

class_vector =  zeros(length(idx_vector),1);

for i_prof = 1:1:length(idx_vector)    
    % peaks = struct2array(spks(idx_vector(i_prof))); % will use direc biin matrix, line below
    peaks = find(global_events_bin(idx_vector(i_prof),:));
    if not(isempty(peaks))        
        n_pks = length(peaks);        
        % [prom_l,prom_r,prom_max,wact,spk_u,dstim,class] = getProminence(peaks(:,1),profiles,tstim,0);
        if ((n_pks>1)&&(global_events(idx_vector(i_prof),peaks(1))<global_events(idx_vector(i_prof),peaks(2)))&&(peaks(2)<=(floor(tstim/Ts)+th_vals.th1))&&en_art_3) %New criteria
            prom_vector(i_prof) = global_events(idx_vector(i_prof),peaks(2));
            act_vector(i_prof) = peaks(2)*Ts;                    
        else
            prom_vector(i_prof) = global_events(idx_vector(i_prof),peaks(1));
            act_vector(i_prof) = peaks(1)*Ts;        
            % prom_vector(i_prof) = prom_l(1);
            % act_vector(i_prof) = spk_u(1,3);
        end
    else
        prom_vector(i_prof) = 0;
        act_vector(i_prof) = 0;
    end
end

prom_max = max(prom_vector);

th_prom_percentage = 0.15;
th_prom = th_prom_percentage*prom_max;
t=Ts.*[0:1:n_frame-1];
class_vector = prom_vector >= th_prom;

%% PKT NCLS class scheme
pclass_u = unique(class_vector);
n_classes = length(pclass_u);
string_class = cell(1,n_classes);

for i_ncls4=1:1:length(pclass_u)
    NCLS(i_ncls4).cVector = idx_vector(find(class_vector==pclass_u(i_ncls4))');
    disp(strcat("size(Class ",num2str(pclass_u(i_ncls4))," ) = ",num2str(size(NCLS(i_ncls4).cVector))))
    string_class(i_ncls4) = {strcat("Class2.0. ",num2str(pclass_u(i_ncls4)),": ",values(cdict,{pclass_u(i_ncls4)}))};
    % NCLS(i_ncls4).label = strcat("SubClass ",num2str(pclass_u(i_ncls4),'%02u'));
    NCLS(i_ncls4).label = string_class(i_ncls4);
    NCLS(i_ncls4).color = 255.*colors(pclass_u(i_ncls4)+1,1:3);
    NCLS(i_ncls4).title = strcat("Class2.0.",num2str(pclass_u(i_ncls4)),": ",values(cdict,{pclass_u(i_ncls4)}));
    NCLS(i_ncls4).prof_type = mean(profiles([NCLS(i_ncls4).cVector],:),1);
end

NCSL2_0= NCLS;

% %% Scanning different threshold values for selection based on expert
% criteria BEGIN
% 
th_prom_vector = 0.1; %[0.1,0.15,0.2,0.25,0.30,0.35,0.40,0.50];
respath = "C:\Users\marci\OneDrive\Documentos\MATLAB\Results";


% for i_prom = 1:1:length(th_prom_vector)
% 
%     th_prom_percentage = th_prom_vector(i_prom);
%     th_prom = th_prom_percentage*prom_max;
%     t=Ts.*[0:1:n_frame-1];
%     prom_vector_bin1 = prom_vector < th_prom;
%     prom_vector_bin2 = prom_vector >= th_prom;
% 
%     fclass=figure;
%     subplot(2,1,1)
%     if sum(prom_vector_bin1)>0
%         plot(t,profiles(idx_vector(prom_vector_bin1),:),'r')
%         title("Class2.0.0 ")
%         xlim([t(1) t(end)])
%         ylim([0.5 1.5])
%         xline(tstim,'-',{strcat('Tstim=',num2str(tstim))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%     end
%     subplot(2,1,2)
%     if sum(prom_vector_bin2)>0
%         plot(t,profiles(idx_vector(prom_vector_bin2),:),'b')
%         title("Class2.0.1 ")
%         xlim([t(1) t(end)])
%         ylim([0.5 1.5])
%         xline(tstim,'-',{strcat('Tstim=',num2str(tstim))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');       
%     end
%     sgtitle(strcat("Class2-Subclass0--th_prom = ",num2str(th_prom_percentage)))
%     saveas(fclass,fullfile(respath,strcat("IMG_A_",num2str(th_prom_percentage*100,'%3d'),"_",num2str(i_prom))),'png')
% 
%     if sum(prom_vector_bin1)>0
%         idx2plot = idx_vector(prom_vector_bin1);
%         for i_c0=1:1:length(idx2plot)
%             i_c0
%             fig_c0 = figure;        
%             plot(t,profiles(idx2plot(i_c0),:),'r')
%             hold on
%             if act_vector(idx_vector==idx2plot(i_c0)) ~= 0
%                 % plot(t(act_vector(idx2plot(i_c0))),profiles(idx2plot(i_c0),act_vector(idx2plot(i_c0))),'r*')
%                 plot(act_vector(find(idx_vector==idx2plot(i_c0))),profiles(idx2plot(i_c0),floor(act_vector(find(idx_vector==idx2plot(i_c0)))/Ts)),'r*')
%             end
%             xlim([t(1) t(end)])
%             ylim([0.5 1.5])
%             title(strcat("Class 2.0.0 - Th=",num2str(th_prom_percentage),"Prof=",num2str(idx2plot(i_c0))))
%             xline(tstim,'-',{strcat('Tstim=',num2str(tstim))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%             xline(tstim+Ts*(th_vals.th1),'-',{strcat('th1=Tstim+',num2str(Ts*(th_vals.th1)))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%             xline(tstim+Ts*(th_vals.th1+th_vals.th2),'-',{strcat('th2=Tstim+',num2str(Ts*(th_vals.th1+th_vals.th2)))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%             xline(tstim+Ts*(th_vals.th1+th_vals.th2+th_vals.th3),'-',{strcat('th3=Tstim+',num2str(Ts*(th_vals.th1+th_vals.th2+th_vals.th3)))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%             yline(mean(profiles(idx2plot(i_c0),1:floor(tstim/Ts))),'--')
%             hold off
%             saveas(fig_c0,fullfile(respath,strcat("IMG_B_",num2str(th_prom_percentage*100,'%3d'),"_",num2str(i_c0))),'png')
%         end
%     end
% 
%     if sum(prom_vector_bin2)>0
%         idx2plot = idx_vector(prom_vector_bin2);    
%         for i_c1=1:1:length(idx2plot)
%             fig_c1 = figure;        
%             plot(t,profiles(idx2plot(i_c1),:),'b')
%             hold on
%             if act_vector(idx_vector==idx2plot(i_c1)) ~= 0
%                 plot(act_vector(find(idx_vector==idx2plot(i_c1))),profiles(idx2plot(i_c1),floor(act_vector(find(idx_vector==idx2plot(i_c1)))/Ts)),'b*')
%                 % plot(t(act_vector(idx2plot(i_c1))),profiles(idx2plot(i_c1),floor(act_vector(idx2plot(i_c1))/Ts)),'b*')
%             end
%             xlim([t(1) t(end)])
%             ylim([0.5 1.5])        
%             title(strcat("Class 2.0.1 - Th=",num2str(th_prom_percentage),"Prof=",num2str(idx2plot(i_c1))))
%             xline(tstim,'-',{strcat('Tstim=',num2str(tstim))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%             xline(tstim+Ts*(th_vals.th1),'-',{strcat('th1=Tstim+',num2str(Ts*(th_vals.th1)))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%             xline(tstim+Ts*(th_vals.th1+th_vals.th2),'-',{strcat('th2=Tstim+',num2str(Ts*(th_vals.th1+th_vals.th2)))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%             xline(tstim+Ts*(th_vals.th1+th_vals.th2+th_vals.th3),'-',{strcat('th3=Tstim+',num2str(Ts*(th_vals.th1+th_vals.th2+th_vals.th3)))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
%             hold off
%             saveas(fig_c1,fullfile(respath,strcat("IMG_C_",num2str(th_prom_percentage*100,'%3d'),"_",num2str(i_c1,'%03d'))),'png')
%         end
%     end
% 
%     genReport(respath,strcat(exp_name,"_th_",num2str(th_prom_percentage*100,'%3d')),[],[])    
% 
%     filePattern = fullfile(respath, '*.png'); % Change to whatever pattern you need.
%     theFiles = dir(filePattern);
%     for k = 1 : length(theFiles)
%         baseFileName = theFiles(k).name;
%         fullFileName = fullfile(respath, baseFileName);
%         fprintf(1, 'Now deleting %s\n', fullFileName);
%         delete(fullFileName);
%     end
% end
%% Scanning different threshold values for selection based on expert
% criteria END




