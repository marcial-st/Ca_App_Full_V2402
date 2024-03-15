function [CELLS]=cellsCleanOverlap(CELLS)
% tlabel=[CELLS.track_label]; It delivers only uint8, useless
% tlabel=double(zeros(length(CELLS),1));
% for kk=1:1:length(CELLS)
%     tlabel(kk)=CELLS(kk).track_label;
% end
% tlabel_u=unique(tlabel);
% tlabel_frec = histcounts(tlabel,tlabel_u);
% [~,idx_frec]=find(tlabel_frec~=1);
% for k1=1:1:length(idx_frec)
%     if ((tlabel_u(idx_frec(k1)~=0)||(idx_frec(k1)~=9999))
%     [~,idx_tl]=find(tlabel==tlabel_u(idx_frec(k1)));
%     for k2=1:1:length(idx_tl)
%         CELLS(idx_tl(k2)).status = 'OFFLINE';
%     end
%     end
for kk=1:1:length(CELLS)        
    for jj=1:1:length(CELLS)
        if (jj~=kk)
            if (CELLS(kk).bbox(1:2) == CELLS(jj).bbox(1:2))
                CELLS(kk).status = 'OFFLINE';
                CELLS(jj).status = 'OFFLINE';
            end
        end
    end
end
