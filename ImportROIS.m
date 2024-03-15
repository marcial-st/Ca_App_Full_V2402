T = readtable('C:\Users\MST\Google Drive\CaExp_2020\191017e1\ROIS.xlsx');
pos = table2array(T);
[n_pos,~] = size(pos);

load imeb
imebl = logical(imeb);
[label,n_label] = bwlabel(imebl);
stats = regionprops(imebl,'Centroid');
centroids = reshape(floor([stats.Centroid]),2,n_label)';

imres = insertMarker(imeb,pos,'color','red','size',10);
figure;imshow(imres)

for i = 1:1:n_pos
    dif = pos(1,:).*ones(n_label,2)-centroids;
    dist = diag(sqrt(dif*dif'));
    [min_val,min_idx] = min(dist);
%     pos_update = centroids(min_idx,:);
    pos_update = centroids(min_idx,:);
end

