function [prom_l,prom_r,prom_max,wact,spk_u] = getFeatures(spk_u,ca_prof_smooth)

minimos = find(islocalmin(ca_prof_smooth));
actividades = floor(spk_u(:,1)/3);

prom_l = zeros(size(actividades));
prom_r = prom_l;
prom_max = prom_l;
wact = prom_l;

for k=1:1:length(actividades)
    
    if ismember(actividades(k),minimos)        
        aux = sort(minimos);
    else        
        aux = sort([minimos,actividades(k)]);
    end
    
    widx = find(aux == actividades(k));
    
    if ((widx ~= 1)&&(widx ~= length(aux)))
        epsilon = floor(abs(aux(widx-1)-aux(widx)));        
        prom_l(k) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon))-ca_prof_smooth(aux(widx-1));
        prom_r(k) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon))-ca_prof_smooth(aux(widx+1));
        prom_max(k) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon))-max([ca_prof_smooth(aux(widx+1)),ca_prof_smooth(aux(widx-1))]);
        wact(k) = aux(widx+1)-aux(widx-1);
        spk_u(k,3) = find(ca_prof_smooth== max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon)))*3;
        spk_u(k,4) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon));
%         prom_l(k) = max(ca_prof_smooth(aux(widx-1):aux(widx+1)))-ca_prof_smooth(aux(widx-1));
%         prom_r(k) = max(ca_prof_smooth(aux(widx-1):aux(widx+1)))-ca_prof_smooth(aux(widx+1));
%         prom_max(k) = max(ca_prof_smooth(aux(widx-1):aux(widx+1)))-max([ca_prof_smooth(aux(widx+1)),ca_prof_smooth(aux(widx-1))]);
%         wact(k) = aux(widx+1)-aux(widx-1);
%         spk_u(k,3) = find(ca_prof_smooth== max(ca_prof_smooth(aux(widx-1):aux(widx+1))))*3;
%         spk_u(k,4) = max(ca_prof_smooth(aux(widx-1):aux(widx+1)));
    end                                                                 
end