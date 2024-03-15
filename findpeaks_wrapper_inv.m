 function [pks_all,pk_max_loc_f,pks_w_all,pks_prom_all,tri_vector] = findpeaks_wrapper_inv(profile_smooth,visual_en,i_prof)
Ts = 3;
min_sep = 0;
% profile_smooth = profiles_smooth(90,:);
tri_vector = zeros(1,length(profile_smooth));
t=(0:1:length(profile_smooth)-1)*Ts;
s = get(0, 'ScreenSize');
[pk_max_loc_f,pks_all] = islocalmax(profile_smooth,'MinSeparation',min_sep);
pks_all = profile_smooth(pk_max_loc_f);
[pk_min_loc,~] = islocalmin(profile_smooth,'MinSeparation',min_sep);

pk_max_loc_f = find(pk_max_loc_f);
pk_min_loc = find(pk_min_loc);

pks_prom_all = zeros(length(pk_max_loc_f),1);
pks_w_all = zeros(length(pk_max_loc_f),1);

if visual_en
    figure('Position',[0 0 s(3) s(4)]);plot(t,profile_smooth);hold on
end

for k=1:1:length(pk_max_loc_f)
    aux_min = find(pk_min_loc>pk_max_loc_f(k));
    if ~isempty(aux_min)
        pks_prom_all(k) = profile_smooth(pk_max_loc_f(k))-profile_smooth(pk_min_loc(aux_min(1)));
        tri_vector(pk_max_loc_f(k)) = 1;
        tri_vector(pk_min_loc(aux_min(1))) = -1;
        pks_w_all(k) = (pk_min_loc(aux_min(1))-pk_max_loc_f(k))*Ts;
        if visual_en
            colorm = (0.85-0.25).*rand(1,3) + 0.25;
            plot((pk_max_loc_f(k)-1)*Ts,profile_smooth(pk_max_loc_f(k)),'Marker','p','MarkerEdgeColor',colorm,'MarkerSize',10)
            plot((pk_min_loc(aux_min(1))-1)*Ts,profile_smooth(pk_min_loc(aux_min(1))),'Marker','^','MarkerEdgeColor',colorm,'MarkerSize',5)
            text((pk_max_loc_f(k)-1)*Ts,profile_smooth(pk_max_loc_f(k))+0.01,num2str(pks_prom_all(k),'%0.5f'),'Color',colorm,'FontSize',10,'Rotation',90);        
        end
    elseif(k == length(pk_max_loc_f))
        if profile_smooth(pk_max_loc_f(k))>profile_smooth(end)
            pks_prom_all(k) = profile_smooth(pk_max_loc_f(k))-profile_smooth(end);
            pks_w_all(k) = (length(profile_smooth)-pk_max_loc_f(k))*Ts;
            tri_vector(pk_max_loc_f(k)) = 1;
            tri_vector(end) = -1;
            if visual_en
                colorm = (0.85-0.25).*rand(1,3) + 0.25;
                plot((pk_max_loc_f(k)-1)*Ts,profile_smooth(pk_max_loc_f(k)),'Marker','p','MarkerEdgeColor',colorm,'MarkerSize',10)
                plot(t(end),profile_smooth(end),'Marker','^','MarkerEdgeColor',colorm,'MarkerSize',5)
                text((pk_max_loc_f(k)-1)*Ts,profile_smooth(pk_max_loc_f(k))+0.01,num2str(pks_prom_all(k),'%0.5f'),'Color',colorm,'FontSize',10,'Rotation',90);
            end
        end
    else
        disp(strcat("************** Not detected Prominence pk idx=",num2str(pk_max_loc_f(k))));
    end
end

% if pk_min_loc(end)>= pk_max_loc_f(end)
%     if profile_smooth(end)>profile_smooth(pk_min_loc(end))
%         pks_prom_all(k+1) = profile_smooth(end)-profile_smooth(pk_min_loc(end));
%         pks_w_all(k+1) = (length(profile_smooth)-pk_min_loc(end))*Ts;
%         tri_vector(end) = 1;
%         tri_vector(pk_min_loc(end)) = -1;
%         if visual_en
%             plot(t(end),profile_smooth(end),'Marker','o','MarkerEdgeColor','r')
%         end
%         pks_all(k+1) = profile_smooth(end);
%         pk_max_loc_f(k+1) = length(profile_smooth);
%     end
% else
%     if profile_smooth(end)>profile_smooth(pk_max_loc_f(end))
%         pks_prom_all(k+1) = profile_smooth(end)-profile_smooth(pk_max_loc_f(end));
%         pks_w_all(k+1) = (length(profile_smooth)-pk_max_loc_f(end))*Ts;        
%         tri_vector(end) = 1;
%         tri_vector(pk_max_loc_f(end)) = -1;
%         if visual_en
%             plot(t(end),profile_smooth(end),'Marker','o','MarkerEdgeColor','r')
%         end
%         pks_all(k+1) = profile_smooth(end);
%         pk_max_loc_f(k+1) = length(profile_smooth);
%     end
% end


if visual_en
    ylim([0.7 1.5]);ylabel("ratio",'FontSize',14)
    xlim([0 Ts*length(profile_smooth)]);xlabel("seconds",'FontSize',14)
    title(strcat("Profile: ",num2str(i_prof)),'FontSize',22)    
end