function [sim_cell]=syntheticClassAssign(sim_cell,inj_x,color_matrix,probabilities_x,sim_col,sim_row,m_canv,n_canv,sim_cw,sim_ch,pix2um)
simc_numcells = max(size(sim_cell));
%  Fuzzy approach begin
% fis=readfis('fuzzy_classifier_v200416.fis');
%  Fuzzy approach end

factor_cell_canvas_x=n_canv/sim_col;

for i=1:1:simc_numcells            
    dist2inj = abs(sim_cell(i).params(1)-inj_x)/factor_cell_canvas_x;
    sim_cell(i).dist2inj = dist2inj*sim_cw*pix2um;
    selcol = floor(dist2inj)    
    if selcol == 0 
        selcol=1;
    end
    if selcol>length(probabilities_x)
        selcol=length(probabilities_x);
    end
    classi = randsample([1 2 3 4],1,true,probabilities_x(:,selcol)');
    %  Fuzzy approach begin
    % input = abs(sim_cell(i).params(1)-inj_x);
    % classi = evalfis(fis,input);  % Fuzzy approach
    % if classi >= 0.5
    %     sim_cell(i).class = 2;        
    % else
    %     sim_cell(i).class = 1;
    % end
    sim_cell(i).class = classi;
    sim_cell(i).color = color_matrix(sim_cell(i).class,:);
end