function sim_cell = scaleRandomClasses(sim_cell);

n_scells = length(sim_cell);

aux_icp_mean_vector = [];

kv = 1;
for k=1:1:n_scells
    if (sim_cell(k).class == 3 )
        aux_icp_mean_vector(kv) = mean(sim_cell(k).icp);
        kv = kv +1;
    end
end
m_icps = mean(aux_icp_mean_vector); 

for k=1:1:n_scells
    if ((sim_cell(k).class == 4 ))
        sim_cell(k).icp = sim_cell(k).icp.*0.2+1.2*m_icps;
    end
    if ((sim_cell(k).class == 2 ))
        sim_cell(k).icp = sim_cell(k).icp.*0.2+1.1*m_icps;
    end    
end
