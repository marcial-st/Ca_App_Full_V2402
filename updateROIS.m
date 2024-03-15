function [pos_updated] = updateROIS(pos,LBase)

[n_pos,~] = size(pos);
labels = unique(LBase);
n_label = length(labels)-1;
imebl = LBase > 0;

stats = regionprops(imebl,'Centroid');
centroids = reshape(floor([stats.Centroid]),2,n_label)';

dist_min = zeros(n_pos,2);
for i = 1:1:n_pos
    dif = pos(i,:).*ones(n_label,2)-centroids;
    dist = diag(sqrt(dif*dif'));
    [min_val,min_idx] = min(dist);    
    
    dist_min(i,1) = min_val;
    dist_min(i,2) = min_idx;         
end

dist_th = mean(dist_min(:,1))+3*std(dist_min(:,1));
pos_updated = zeros(size(pos));
label_updated = zeros(n_pos,1);

for i = 1:1:n_pos
    if dist_min(i,1) <= dist_th
        pos_updated(i,:) = centroids(dist_min(i,2),:);       
        label_updated(i) = LBase(centroids(dist_min(i,2),2),centroids(dist_min(i,2),1));
    else
        pos_updated(i,:) = [-1,-1];
        label_updated(i) = 0;
    end
end

i_upd = 1;
i_con = 1;
for i = 1:1:n_pos
    aux= find(label_updated == label_updated(i));
    if length(aux)>1
        pos_con(i_con,:) = pos(i,:);
        pos_updated(i,:) = [-1,-1];
        i_con = i_con+1;
    else
        pos_upd(i_upd,:) = pos_updated(i,:);        
        i_upd = i_upd+1;
    end
end
[~,n_lbl] = size(LBase);
rois_original = figure;
axes(rois_original);
imres = insertMarker(double(imebl),pos,'x','color','green','size',5);
imres = insertMarker(imres,pos_upd,'+','color','blue','size',5);    
imres = insertMarker(imres,pos_upd,'o','color','blue','size',5);    
imres = insertMarker(imres,pos_con,'x','color','red','size',5);    
imres = insertMarker(imres,pos_con,'o','color','red','size',5);    
label_str = {'Original position','Updated position','Conflict'};
position = [n_lbl-120 25 0 0;n_lbl-120 50 0 0;n_lbl-120 75 0 0];
imres = insertObjectAnnotation(imres,'rectangle',position,label_str,'TextBoxOpacity',0.8,'FontSize',12,'Color',{'green','blue','red'});
imshow(imres);
title('ROIs location')
