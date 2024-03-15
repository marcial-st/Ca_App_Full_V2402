function sytheticFramePlot(sim_cell,sim_canv)
% figure;
imshow(sim_canv);
hold on;
n_scells = max(size(sim_cell));
for k=1:1:n_scells        
        rectangle('Position',sim_cell(k).params,'Curvature',[1 1],'EdgeColor',sim_cell(k).color,"FaceColor",sim_cell(k).color) % Draw s-cell> per cell   
end

hold off;