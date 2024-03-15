function [sim_cell,sim_canv] = generateSyntheticFrame(sim_row,sim_col)
% 1) function [sim_cell]=generateSyntheticFrame(sim_row,sim_col)
% donde  
%   sim_cell : estructura de datos de las células sintéticas   
%    sim_row : número de renglones de células sintéticas
%    sim_col : número de columnas de células sintéticas
% 
% 2) function [synthetic_image]=plotSyntheticFrame(sim_cell)  
% synthetic_image: imagen/frame sintética
simc_numcells = sim_col*sim_row;

sim_cw = 40; %% Maximum width of r-cell bounding box based on exp. variability
sim_ch = 20; %% Maximum heigt of r-cell bounding box based on exp. variability

mesh_rows = 0:2*sim_ch:sim_row*2*sim_ch; %% Mesh vectors
mesh_cols = 0:2*sim_cw:sim_col*2*sim_cw; % 2*MaxWidth & 2*MaxHeight to avoid collisions between s-cells

sim_canv_x = sim_col*2*sim_cw;
sim_canv_y = sim_row*2*sim_ch;
sim_canv = uint8(zeros(sim_canv_y,sim_canv_x)); % empty canvas generation
sim_canv(:,:) = 100;                                        % gray background just for visualization
sim_canv(mesh_rows(2:end),:) = 255;                         %% Mesh lines
sim_canv(:,mesh_cols(2:end)) = 255;                         %% Mesh bars

simc_x = round(sim_cw.*rand(simc_numcells,1));              %% Random position vectors x&y of the i-th s-cell in its own cell (Uniform distribution)
simc_y = round(sim_ch.*rand(simc_numcells,1));              
% TODO: code a function to extract the parameters automatically regardless
% the type of tissue-cell
% simc_w = simc_std_w.*randn(simc_numcells,1)+simc_mean_w;   %% Random geometry vectors w&h of the i-th s-cell (Normal distribution)
simc_w = 4.4255.*randn(simc_numcells,1)+19.0270;
simc_h = simc_w./(0.3072*randn(simc_numcells,1)+1.5060);   

% figure;
% imshow(sim_canv);
% hold on;
k=0; % s-cells vectors parameters index

for i=1:1:sim_row
    for j=1:1:sim_col                
        k = k+1;
        sim_cell(k).params = [mesh_cols(j)+simc_x(k) mesh_rows(i)+simc_y(k) simc_w(k) simc_h(k)]; %% Build structure for each s-cell
%         rectangle('Position',sim_cell(k).params,'Curvature',[1 1],'EdgeColor','r',"FaceColor",'r') % Draw s-cell> per cell
    end
end

% hold off;