function [sim_image] = sytheticFramPlot(sim_cell,sim_canv)
% figure;
imshow(sim_canv);
hold on;
k=0; % s-cells vectors parameters index

for i=1:1:sim_row
    for j=1:1:sim_col                
        k = k+1;        
        rectangle('Position',sim_cell(k).params,'Curvature',[1 1],'EdgeColor','r',"FaceColor",'r') % Draw s-cell> per cell
    end
end

hold off;