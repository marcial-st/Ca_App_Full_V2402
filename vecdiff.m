function [difv] = vecdiff(vec1,vec2)
lv1=length(vec1);
lv2=length(vec2);
dindx = 1;
if lv1==lv2
    for i=1:1:lv1
        if vec1(i) ~= vec2(i)
            difv(dindx) = i;
            dindx = dindx+1;
        end
    end
else
    warndlg('Vectors must have the same length')
end