function [global_events_bin_updated,global_events_updated] = UpdateGlobalEvents_th(profiles_smooth,global_events_bin_vector,global_events_vector,global_events_w_vector,Ts,th_p)
global_events_bin_updated = zeros(size(global_events_bin_vector));
global_events_updated = zeros(size(global_events_vector));
idx_e = find(global_events_bin_vector);
if ~isempty(idx_e)
    for i=1:1:length(idx_e)
        if global_events_w_vector(idx_e(i))~=0                        
            idx_new = idx_e(i)-floor(global_events_w_vector(idx_e(i))/Ts);
            idx_max = idx_e(i);
            idx_min = idx_e(i)-floor(global_events_w_vector(idx_e(i))/Ts);
            if idx_min>=0
                if idx_min==0
                    idx_min = 1;
                end
                
                prom = profiles_smooth(idx_max)-profiles_smooth(idx_min);                
                valth = profiles_smooth(idx_min)+prom*th_p;                
                profile_valth = profiles_smooth(idx_min:idx_max)-valth;
                idx_new = find(diff(profile_valth>=0),1)+idx_min;                                 
                
                global_events_bin_updated(idx_e(i)) = 0;
                global_events_bin_updated(idx_new) = 1;
                
                global_events_updated(idx_new) = global_events_vector(idx_e(i));
                global_events_vector(idx_e(i)) = 0;
                
%                 figure;
%                 plot(profiles_smooth,'bo-');hold on;
%                 plot(idx_max,profiles_smooth(idx_max),'r*','MarkerSize',12)
%                 plot(idx_min,profiles_smooth(idx_min),'r*','MarkerSize',12)
%                 yline(profiles_smooth(idx_min),'r--')
%                 plot(idx_new,profiles_smooth(idx_new),'mo','MarkerSize',14)
%                 plot(idx_new,profiles_smooth(idx_new),'m+','MarkerSize',14)
%                 plot([idx_new idx_new],[profiles_smooth(idx_min) profiles_smooth(idx_new)],'m:')
%                 plot([1 idx_new],[profiles_smooth(idx_new) profiles_smooth(idx_new)],'m:')
%                 title(strcat(num2str(100*th_p),"%"))
                
            else
                disp("UpdateGlobalEvents Warning: idx<=0, not able to update this location")
            end
        end
    end
end