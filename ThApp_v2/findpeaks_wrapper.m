 function [pks_all,pk_max_loc_f,pks_w_all,pks_prom_all] = findpeaks_wrapper(profile_smooth)
Ts = 3;
% profile_smooth = profiles_smooth(90,:);
% t=(0:1:length(profile_smooth)-1)*Ts;

[pk_max_loc_f,pks_all] = islocalmax(profile_smooth);
pks_all = profile_smooth(pk_max_loc_f);
[pk_min_loc,~] = islocalmin(profile_smooth);

pk_max_loc_f = find(pk_max_loc_f);
pk_min_loc = find(pk_min_loc);

pks_prom_all = zeros(length(pk_max_loc_f),1);
pks_w_all = zeros(length(pk_max_loc_f),1);

for k=1:1:length(pk_max_loc_f)
    aux_min = find(pk_min_loc<pk_max_loc_f(k));
    if ~isempty(aux_min)
        pks_prom_all(k) = profile_smooth(pk_max_loc_f(k))-profile_smooth(pk_min_loc(aux_min(end)));
        pks_w_all(k) = (pk_max_loc_f(k)-pk_min_loc(aux_min(end)))*Ts;
    elseif(k == 1)
        pks_prom_all(k) = profile_smooth(pk_max_loc_f(k))-profile_smooth(1);
        pks_w_all(k) = (pk_max_loc_f(k))*Ts;        
    else
        disp(strcat("************** Not detected Prominence pk idx=",num2str(pk_max_loc_f(k))));
    end
end