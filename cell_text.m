%% Thesis
function [text_m]=cell_text(CELLS)
[r,c]=size(CELLS);
text_m=cell(c,1);
for i=1:1:c
    text_m(i)=strcat({'Cell '},CELLS(i).user_label,{':  '},CELLS(i).status);
end