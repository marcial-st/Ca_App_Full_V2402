function [NCLS16,NCLS4,th_vals] = subclass2_classifier(spks,idx_vector,Ts,tstim,profiles,global_events,global_events_bin)

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

th_vals.th1 = 20;
th_vals.th2 = 10;
th_vals.th3 = 10;

th1=th_vals.th1*Ts;
th2=th_vals.th2*Ts;
th3=th_vals.th3*Ts;

cindex4 = [0,1,2,3];
clabels4 = ["th1"; %0
            "th2"; %1
            "th3"; %2
            ">th3"];%3
    
cdict4 = containers.Map(cindex4,clabels4);


cindex16 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
clabels16 = [">th3=0,th3=0,th2=0,th1=0"; %0
             ">th3=0,th3=0,th2=0,th1=1"; %1
             ">th3=0,th3=0,th2=1,th1=0"; %2
             ">th3=0,th3=0,th2=1,th1=1"; %3
             ">th3=0,th3=1,th2=0,th1=0"; %4
             ">th3=0,th3=1,th2=0,th1=1"; %5
             ">th3=0,th3=1,th2=1,th1=0"; %6
             ">th3=0,th3=1,th2=1,th1=1"; %7
             ">th3=1,th3=0,th2=0,th1=0"; %8
             ">th3=1,th3=0,th2=0,th1=1"; %9
             ">th3=1,th3=0,th2=1,th1=0"; %10
             ">th3=1,th3=0,th2=1,th1=1"; %11
             ">th3=1,th3=1,th2=0,th1=0"; %12
             ">th3=1,th3=1,th2=0,th1=1"; %13
             ">th3=1,th3=1,th2=1,th1=0"; %14
             ">th3=1,th3=1,th2=1,th1=1"];%15
cdict16 = containers.Map(cindex16,clabels16);

[n_prof,n_frame]=size(profiles);

class_vector =  zeros(length(idx_vector),1);

class_matrix = zeros(length(idx_vector),4);
class_matrix_cdec = zeros(length(idx_vector),1);

for i_prof = 1:1:length(idx_vector)
    % peaks = struct2array(spks(idx_vector(i_prof))); % will use direc biin matrix, line below
    peaks = find(global_events_bin(idx_vector(i_prof),:));
    if not(isempty(peaks))
        peaks = Ts*peaks;
        [n_pks,~] = size(peaks);        

        if (peaks(1,1)>=(tstim) && peaks(1,1)<(tstim+th1))
            class_vector(i_prof) = 0;           
        elseif (peaks(1,1)>=(tstim+th1) && peaks(1,1)<(tstim+th1+th2))
            class_vector(i_prof) = 1;
        elseif (peaks(1,1)>=(tstim+th1+th2) && peaks(1,1)<(tstim+th1+th2+th3))
            class_vector(i_prof) = 2;
        elseif(peaks(1,1)>=(tstim+th1+th2+th3))
            class_vector(i_prof) = 3;
        end
        
        for i_pk=1:1:n_pks
            if (peaks(i_pk,1)>=(tstim) && peaks(i_pk,1)<(tstim+th1))
                class_matrix(i_prof,1)=1;                
            elseif (peaks(i_pk,1)>=(tstim+th1) && peaks(i_pk,1)<(tstim+th1+th2))
                class_matrix(i_prof,2)=1;
            elseif (peaks(i_pk,1)>=(tstim+th1+th2) && peaks(i_pk,1)<(tstim+th1+th2+th3))
                class_matrix(i_prof,3)=1;
            elseif (peaks(i_pk,1)>=(tstim+th1+th2+th3))
                class_matrix(i_prof,4)=1;
            end
        end
    end          
end
for i_matrix = 1:1:length(idx_vector)
    class_matrix_dec(i_matrix) = bin2dec(num2str(class_matrix(i_matrix,:)));
end

%% PKT NCLS 4 class scheme
pclass_u = unique(class_vector);
n_classes = length(pclass_u);
string_class = cell(1,n_classes);

for i_ncls4=1:1:length(pclass_u)
    NCLS(i_ncls4).cVector = idx_vector(find(class_vector==pclass_u(i_ncls4))');
    disp(strcat("size(Class2S ",num2str(pclass_u(i_ncls4))," ) = ",num2str(size(NCLS(i_ncls4).cVector))))
    string_class(i_ncls4) = {strcat("Class2. ",num2str(pclass_u(i_ncls4)),": ",values(cdict4,{pclass_u(i_ncls4)}))};
    % NCLS(i_ncls4).label = strcat("SubClass ",num2str(pclass_u(i_ncls4),'%02u'));
    NCLS(i_ncls4).label = string_class(i_ncls4);
    NCLS(i_ncls4).color = 255.*colors(pclass_u(i_ncls4)+1,1:3);
    NCLS(i_ncls4).title = strcat("Class2: Subclass ",num2str(pclass_u(i_ncls4)),": ",values(cdict4,{pclass_u(i_ncls4)}));
    NCLS(i_ncls4).prof_type = mean(profiles([NCLS(i_ncls4).cVector],:),1);
end 
NCLS4 = NCLS;

%% PKT NCLS 16 class scheme
pclass_u = unique(class_matrix_dec);
n_classes = length(pclass_u);
string_class = cell(1,n_classes);

for i_ncls4=1:1:length(pclass_u)
    NCLS(i_ncls4).cVector = idx_vector(find(class_matrix_dec==pclass_u(i_ncls4))');
    disp(strcat("size(Class ",num2str(pclass_u(i_ncls4))," ) = ",num2str(size(NCLS(i_ncls4).cVector))))
    string_class(i_ncls4) = {strcat("Class ",num2str(pclass_u(i_ncls4)),": ",values(cdict16,{pclass_u(i_ncls4)}))};
    % NCLS(i_ncls4).label = strcat("SubClass ",num2str(pclass_u(i_ncls4),'%02u'));
    NCLS(i_ncls4).label = string_class(i_ncls4);
    NCLS(i_ncls4).color = 255.*colors(pclass_u(i_ncls4)+1,1:3);
    NCLS(i_ncls4).title = strcat("Class2: Subclass ",num2str(pclass_u(i_ncls4)),": ",values(cdict16,{pclass_u(i_ncls4)}));
    NCLS(i_ncls4).prof_type = mean(profiles([NCLS(i_ncls4).cVector],:),1);
end 
NCLS16 = NCLS;

offset = 0.5;
t=Ts.*[0:1:n_frame-1];

figure;
hold on
for i_c=1:1:length(NCLS4)    
    plot(t,profiles([NCLS4(i_c).cVector],:)+((i_c-1)*offset),'Color',[NCLS4(i_c).color]./255);
    text(t(end)-150,1.2+(i_c-1)*offset,NCLS4(i_c).label,'Color','m','Interpreter','none')
end
xline(tstim,'-',{strcat('Tstim=',num2str(tstim))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
xline(tstim+th1,'-',{strcat('Tstim+',num2str(th1))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
xline(tstim+th1+th2,'-',{strcat('Tstim+',num2str(th1+th2))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
xline(tstim+th1+th2+th3,'-',{strcat('Tstim+',num2str(th1+th2+th3))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
xlim([t(1) t(end)])
xlabel('Time [s]')
title('Class2_Subclasses: 4 scheme','Interpreter','none')
hold off

figure;
hold on
for i_c=1:1:length(NCLS4)    
    plot(t,profiles([NCLS4(i_c).cVector],:)+((i_c-1)*offset),'Color',[200 200 200]./255);
    plot(t,NCLS4(i_c).prof_type+((i_c-1)*offset),'Color',[NCLS4(i_c).color]./255,'LineWidth',2);
    text(t(end)-150,1.2+(i_c-1)*offset,NCLS4(i_c).label,'Color','m','Interpreter','none')
end
xline(tstim,'-',{strcat('Tstim=',num2str(tstim))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
xline(tstim+th1,'-',{strcat('Tstim+',num2str(th1))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
xline(tstim+th1+th2,'-',{strcat('Tstim+',num2str(th1+th2))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
xline(tstim+th1+th2+th3,'-',{strcat('Tstim+',num2str(th1+th2+th3))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
xlim([t(1) t(end)])
xlabel('Time [s]')
title('Class2_Subclasses: 4 scheme','Interpreter','none')
hold off


% figure;
% hold on
% for i_c=1:1:length(NCLS16)    
%     plot(t,profiles([NCLS16(i_c).cVector],:)+((i_c-1)*offset),'Color',[NCLS16(i_c).color]./255);
%     text(t(end)-150,1.2+(i_c-1)*offset,NCLS16(i_c).label,'Color','m','Interpreter','none')
% end
% xline(tstim,'-',{strcat('Tstim=',num2str(tstim))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% xline(tstim+th1,'-',{strcat('Tstim+',num2str(th1))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% xline(tstim+th1+th2,'-',{strcat('Tstim+',num2str(th1+th2))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% xline(tstim+th1+th2+th3,'-',{strcat('Tstim+',num2str(th1+th2+th3))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% xlim([t(1) t(end)])
% xlabel('Time [s]')
% title('Class2_Subclasses: 16 scheme','Interpreter','none')
% hold off
% 
% figure;
% hold on
% for i_c=1:1:length(NCLS16)    
%     plot(t,profiles([NCLS16(i_c).cVector],:)+((i_c-1)*offset),'Color',[200 200 200]./255);
%     plot(t,NCLS16(i_c).prof_type+((i_c-1)*offset),'Color',[NCLS16(i_c).color]./255,'LineWidth',2);
%     text(t(end)-150,1.2+(i_c-1)*offset,NCLS16(i_c).label,'Color','m','Interpreter','none')
% end
% xline(tstim,'-',{strcat('Tstim=',num2str(tstim))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% xline(tstim+th1,'-',{strcat('Tstim+',num2str(th1))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% xline(tstim+th1+th2,'-',{strcat('Tstim+',num2str(th1+th2))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% xline(tstim+th1+th2+th3,'-',{strcat('Tstim+',num2str(th1+th2+th3))},'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% xlim([t(1) t(end)])
% xlabel('Time [s]')
% title('Class2_Subclasses: 16 scheme','Interpreter','none')
% hold off
