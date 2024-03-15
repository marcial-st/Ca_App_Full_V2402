function [prom_l,prom_r,prom_max,wact,spk_u,dstim,class] = getProminence(spk_u,ca_prof_smooth,tstim,na)

minimos = find(islocalmin(ca_prof_smooth));
actividades = floor(spk_u(:,1)/3);

prom_l = zeros(size(actividades));
prom_r = prom_l;
prom_max = prom_l;
wact = prom_l;
dstim = prom_l;
class = prom_l;

for k=1:1:length(actividades)
    
    if ismember(actividades(k),minimos)        
        aux = sort(minimos);
    else        
        aux = sort([minimos,actividades(k)]);
    end
    
    widx = find(aux == actividades(k));
    if na == 0 
        if ((widx ~= 1)&&(widx ~= length(aux)))
            epsilon = floor(abs(aux(widx-1)-aux(widx)));        
            prom_l(k) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon))-ca_prof_smooth(aux(widx-1));
            prom_r(k) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon))-ca_prof_smooth(aux(widx+1));
            prom_max(k) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon))-max([ca_prof_smooth(aux(widx+1)),ca_prof_smooth(aux(widx-1))]);
            wact(k) = 2*epsilon;
            dstim(k) = aux(widx)-tstim;
            spk_u(k,3) = find(ca_prof_smooth== max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon)))*3;
            spk_u(k,4) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon));
            class(k,1) = 'A';
        end                                                                 
    else
        if ((widx ~= 1))
            epsilon = floor(abs(aux(widx-1)-aux(widx)));
            if ((aux(widx)+epsilon) < length(ca_prof_smooth))
                prom_l(k) = max(ca_prof_smooth(aux(widx-1):aux(widx)+epsilon))-ca_prof_smooth(aux(widx-1));
                prom_r(k) = 0;
                prom_max(k) = 0;
                wact(k) = 2*epsilon;
                dstim(k) = aux(widx)-tstim;
                class(k,1) = 'B';
            end
        end
    end
    if ((prom_l(k)==0)&&(wact(k)==0)&&(dstim(k)==0))
        disp("A chis a chis los mariachis")
        na
    end
    
end