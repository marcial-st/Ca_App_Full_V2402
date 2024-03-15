function [sim_cell] = applyCellModelpClass(sim_cell,nNCLS,t)

n_cells = length(sim_cell);
n_class = length(nNCLS);
            
            %"sin_constant"
            % Fss=1; Tss = 1/Fss;
            % n=100*Tss;
            % phi=0;
            % t=0:1/n:2*pi;
            % fluo_max = 0.5;
            % offset = 0.5;
            % icp = offset+fluo_max*sin(Fss*t+phi);
            % sim_cell(i_cell).icp = icp;


for i_cell =1:1:n_cells           
    switch(sim_cell(i_cell).class)
        case 1 %"Class0"
            n_synth_prof = 1;
                        
            offset = (nNCLS(1).fitbasal_params(1)*sim_cell(i_cell).dist2inj)+nNCLS(1).fitbasal_params(2);

            % coeffs_matrix = genCoeffsMatrix(nNCLS(1).fitParams,n_synth_prof); %Commented            
            icp = genProfileFit(nNCLS(1).fitParams,t,nNCLS(1).fitmodel,1,offset);
            fluo_max = 0.5;
                        
            sim_cell(i_cell).icp = icp;    

        case 3 %"Class0"
            coeffs(1) = genProfileFit(nNCLS(3).fitparams_basal,sim_cell(i_cell).dist2inj,nNCLS(3).fitparams_basal_model,0,0);
            % prominence ==> zeta
            coeffs(2) = genProfileFit(nNCLS(3).fitparams_prominence,sim_cell(i_cell).dist2inj,nNCLS(3).fitparams_prominence_model,0,0);
            % t_r ==> w_n
            coeffs(3) = genProfileFit(nNCLS(3).fitparams_peaktime,sim_cell(i_cell).dist2inj,nNCLS(3).fitparams_peaktime_model,0,0);
            % plateau
            coeffs(4) = genProfileFit(nNCLS(3).fitparams_plateau,sim_cell(i_cell).dist2inj,nNCLS(3).fitparams_plateau_model,0,0);
            % Probability of plateau
            coeffs(5) =  nNCLS(3).prob_plateau

            coeffs(6) = nNCLS(3).th_plateu;
            coeffs(7) = nNCLS(3).plateau_range_low    ;    
            coeffs(8) = nNCLS(3).plateau_range_medium ;
            coeffs(9) = nNCLS(3).plateau_range_high   ;
    
            coeffs(10) = nNCLS(3).plateau_prob_low    ;    
            coeffs(11) = nNCLS(3).plateau_prob_medium ;
            coeffs(12) = nNCLS(3).plateau_prob_high   ;

            if (sum(coeffs<0)>0)
                warning(strcat("Got coefficients <0 for d=",num2str(sim_cell(i_cell).dist2inj),", value was limited to 0, please check distance ranges."))
                coeffs
                coeffs(coeffs<0)=0;  % crop negative values
            else
                icp = genProfileFit(coeffs,t,nNCLS(3).fitmodel,0,0);
                % figure;
                % plot(t,icp)
                % xlabel("Time (s)")
                % ylabel("Ratio")
                % title("Class2: Synthetic profile - applyCellModelpClass")
                sim_cell(i_cell).icp = icp;                 
            end                                  

        otherwise % other classes = "random"              
            icp = randn(1,length(t)).*0.1;            
            sim_cell(i_cell).icp = icp;                      
    end

end