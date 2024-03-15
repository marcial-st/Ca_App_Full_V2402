function [sim_cell] = applyCellModel(sim_cell,model_sel,sim_row,sim_col)
n_cells = length(sim_cell);
switch(model_sel)
    case "sin_constant"
        Fss=1; Tss = 1/Fss;
        n=100*Tss;
        phi=0;
        t=0:1/n:2*pi;
        fluo_max = 0.5;
        offset = 0.5;
        icp = offset+fluo_max*sin(Fss*t+phi);
        for i=1:n_cells
            sim_cell(i).icp = icp;
        end
    case "sin_linear_col"
        Fss=1; Tss = 1/Fss;
        n=100*Tss;
        phi=pi/2;
        t=0:1/n:2*pi;
        fluo_max = 0.5;
        offset = 0.5;
        i_row = 1;
        for i=1:n_cells
            phi=(2*pi)*(i_row/sim_row);
            icp = offset+fluo_max*sin(Fss*t+phi);            
            sim_cell(i).icp = icp;
            if i_row == sim_row
                i_row=1;
            else
                i_row = i_row+1;
            end
        end
    case "sin_linear_row"
        Fss=1; Tss = 1/Fss;
        n=100*Tss;
        phi=pi/2;
        t=0:1/n:2*pi;
        fluo_max = 0.5;
        offset = 0.5;
        i_col = 1;
        for i=1:n_cells
            phi=(pi/2)*(i_col/sim_col);
            icp = offset+fluo_max*sin(Fss*t+phi);            
            sim_cell(i).icp = icp;
            if i_col == sim_col
                i_col=1;
            else
                i_col = i_col+1;
            end
        end
    case "random"
          for i=1:n_cells
            icp = rand(1,300);            
            sim_cell(i).icp = icp;
          end
    otherwise
end

