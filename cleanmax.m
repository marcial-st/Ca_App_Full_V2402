function spk_max = cleanmax(spk_max,spk_u);

rep = intersect(spk_max,spk_u);
lrep = length(rep);
if lrep > 0
    for ir =1:1:lrep
        idx = find(spk_max == rep(ir));
        spk_max(idx) = [];
    end
end