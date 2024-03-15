function [sim_cell] = simClassify(sim_cell,inj_x,fis,color_matrix)

n_sim_cells = max(size(sim_cell));
range = fis.Outputs.Range;
n_classes = max(size(fis.Outputs.MembershipFunctions));
thresholds = range(1):range(end)/n_classes:range(end);
thresholds(1)=[];

for k = 1:1:n_sim_cells                    
    input = abs(sim_cell(k).params(1)-inj_x);
    fis_out = evalfis(fis,input);
    
    [~,class] = min(abs(thresholds-fis_out));
    sim_cell(k).class = class;
    sim_cell(k).color = color_matrix(class,:);
%     rectangle('Position',sim_cell(k).params,'Curvature',[1 1],'EdgeColor',color_matrix(sim_cell(k).class,:),"FaceColor",color_matrix(sim_cell(k).class,:)) % Draw s-cell> per cell
end
% rectangle('Position',[inj_x 1 1 sim_row*2*sim_ch-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)