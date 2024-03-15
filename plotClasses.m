function plotClasses(OnCaProfs,NCLS,selclass,t,params,t_sample)
%% Plot the profiles per class
[~,n1]=size(NCLS);
ymax = max(max(OnCaProfs));
ymin = min(min(OnCaProfs));
for i1=1:1:n1
    cVector=NCLS(i1).cVector;
    [~,n2]=size(cVector);
    fclass = figure;
    figure(fclass)
%     titleString=strcat("Class ",num2str(selclass));    
    sgtitle(NCLS(i1).title,'FontName','Times New Roman','Fontsize',16)
%     subplot(1,2,2)
    hold on;
    for i2=1:1:n2
        % plot(t,OnCaProfs(:,cVector(i2))') % Original, all the time range
        fr_start = floor(params.tstart/t_sample);
        fr_end = floor(params.tend/t_sample);
        plot(t(fr_start:fr_end)-t(fr_start),OnCaProfs((fr_start:fr_end),cVector(i2))','LineWidth',1.5) % only analysis time
    end
    if isstruct(params)
        % xline(params.tstart,'--',{'t_{start}'},'LabelVerticalAlignment','bottom','Color',[0.45,0.45,0.45])
        xline(params.tstim-params.tstart,'--',{'t_{inj}'},'LabelVerticalAlignment','top','Color',[0.45,0.45,0.45])
        % xline(params.tend,'--',{'t_{end}'},'LabelVerticalAlignment','bottom','Color',[0.45,0.45,0.45])      
    else
        xline(params(1)*t_sample,'--',{'t_{start}'},'LabelVerticalAlignment','bottom','Color',[0.45,0.45,0.45])
        xline((params(1)+params(2))*t_sample,'--',{'t_{mpeak}'},'LabelVerticalAlignment','bottom','Color',[0.45,0.45,0.45])
        xline(params(5)*t_sample,'--',{'t_{shiftl}'},'LabelVerticalAlignment','bottom','Color',[0.45,0.45,0.45])
        xline(t(end)-params(6)*t_sample,'--',{'t_{shiftr}'},'LabelVerticalAlignment','bottom','Color',[0.45,0.45,0.45])
    end
    set(gca,'Color',[0.85 0.85 0.85])
%     title("Calcium Profiles",'FontName','Times New Roman','Fontsize',10)
    xlabel('Time (seconds)');ylabel('Ratio (F_{340}/F_{380})','Interpreter','tex','FontName','Helvetica');    
    % xlim([t(1) t(end)]); %  Original, all the time range
    xlim([t(fr_start)-t(fr_start) t(fr_end)-t(fr_start)]);
    ylim([ymin ymax]);
    hold off;
%     subplot(1,2,1)        
%     plot(t,NCLS.prof_fit,'Color',[NCLS.color/255],'LineWidth',2)
%     xlabel('Frame');ylabel('Ratio');
%     xlim([0 max(length(t))]);
%     title(strcat("Average Pattern Profile"),'FontName','Times New Roman','Fontsize',10)
end