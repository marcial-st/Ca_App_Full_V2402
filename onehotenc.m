function oht = onehotenc(reflabels,cvector)

oht = zeros(max(reflabels)+1,length(cvector));

for i=1:1:length(cvector)
    oht(cvector(i)+1,i)=1;
end