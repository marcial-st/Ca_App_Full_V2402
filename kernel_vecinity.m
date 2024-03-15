function [kernel_sub] = kernel_vecinity(imDim,vecinity,bbox,label_matrix)
kernel_sub = uint16(zeros(imDim));

if (bbox(1)-vecinity)<1
    v_x1 = 1;
else
    v_x1 = floor(bbox(1)-vecinity);
end
if (bbox(2)-vecinity)<1
    v_y1 = 1;
else
    v_y1 = floor(bbox(2)-vecinity);
end

if (bbox(1)+bbox(3)+vecinity)>imDim(2)
    v_x2 = imDim(2);
else
    v_x2 = floor(bbox(1)+bbox(3)+vecinity);
end
if (bbox(2)+bbox(4)+vecinity)>imDim(1)
    v_y2 = imDim(1);
else
    v_y2 = floor(bbox(2)+bbox(4)+vecinity);
end

kernel_sub(v_y1:v_y2,v_x1:v_x2)=1;
kernel_sub = kernel_sub.*label_matrix;
ks_labels = nonzeros(unique(kernel_sub));
kernel_sub = uint16(zeros(imDim));
for i=1:length(ks_labels)
    kernel_sub = kernel_sub+(ks_labels(i)*uint16((label_matrix == ks_labels(i))));
end