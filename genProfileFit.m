function prof_fit = genProfileFit(coeffs,t,fitmodel,en_offset_override,offset_val)

    switch(fitmodel)
        case 'exp1'
            prof_fit = coeffs(1)*exp(coeffs(2)*t);
        case 'exp2'
            prof_fit = coeffs(1)*exp(coeffs(2)*t) + coeffs(3)*exp(coeffs(4)*t);
        case 'sin3'
            prof_fit = coeffs(1)*sin(coeffs(2)*t+coeffs(3))+coeffs(4)*sin(coeffs(5)*t+coeffs(6))+coeffs(7)*sin(coeffs(8)*t+coeffs(9));
        case 'sin8'
            prof_fit = coeffs(1)*sin(coeffs(2)*t+coeffs(3))+coeffs(4)*sin(coeffs(5)*t+coeffs(6))+coeffs(7)*sin(coeffs(8)*t+coeffs(9))+coeffs(10)*sin(coeffs(11)*t+coeffs(12))+coeffs(13)*sin(coeffs(14)*t+coeffs(15))+coeffs(16)*sin(coeffs(17)*t+coeffs(18))+coeffs(19)*sin(coeffs(20)*t+coeffs(21))+coeffs(22)*sin(coeffs(23)*t+coeffs(24));
        case 'gauss3'
            prof_fit  = coeffs(1).*exp(-((t-coeffs(2))/coeffs(3)).^2) + coeffs(4).*exp(-((t-coeffs(5))/coeffs(6)).^2) + coeffs(7).*exp(-((t-coeffs(8))/coeffs(9)).^2);
        case 'gauss4'
            prof_fit  = coeffs(1).*exp(-((t-coeffs(2))/coeffs(3)).^2) + coeffs(4).*exp(-((t-coeffs(5))/coeffs(6)).^2) + coeffs(7).*exp(-((t-coeffs(8))/coeffs(9)).^2);            
        case 'poly1'
            if en_offset_override
                offset = offset_val;
            else
                offset = coeffs(2);
            end
            prof_fit = coeffs(1)*t+offset;
        case 'poly2'
            if en_offset_override
                offset = offset_val;
            else
                offset = coeffs(3);
            end
            prof_fit = coeffs(1)*t.^2+coeffs(2)*t+offset;
        case 'poly3'    
            if en_offset_override
                offset = offset_val;
            else
                offset = coeffs(4);
            end
            prof_fit = coeffs(1)*t.^3+coeffs(2)*t.^2+coeffs(3)*t+offset;
        case 'poly4'    
            if en_offset_override
                offset = offset_val;
            else
                offset = coeffs(4);
            end
            prof_fit = coeffs(1)*t.^4+coeffs(2)*t.^3+coeffs(3)*t.^2+coeffs(4)*t+offset;            
        case 'poly5'
            if en_offset_override
                offset = offset_val;
            else
                offset = coeffs(6);
            end            
             prof_fit = coeffs(1)*t.^5+coeffs(2)*t.^4+coeffs(3)*t.^3+coeffs(4)*t.^2+coeffs(5)*t+offset;
        case 'secord'
            basal = coeffs(1);
            zeta = abs(log(coeffs(2)))/(sqrt(pi^2+(log(coeffs(2))^2)));           
            w_n  = pi/(coeffs(3)*sqrt(1-zeta^2));
            plateau = coeffs(4);

            % coeffs(6) = NCLS(class_i).th_plateu;
            % coeffs(7) = NCLS(class_i).plateau_range_low    ;    
            % coeffs(8) = NCLS(class_i).plateau_range_medium ;
            % coeffs(9) = NCLS(class_i).plateau_range_high   ;
            % 
            % coeffs(10) = NCLS(class_i).plateau_prob_low    ;    
            % coeffs(11) = NCLS(class_i).plateau_prob_medium ;
            % coeffs(12) = NCLS(class_i).plateau_prob_high   ;            

            en_plateau = randsample([0 1],1,true,[1-coeffs(5) coeffs(5)]);
            disp(strcat("Secord, plateau probability=",num2str(coeffs(5))));

            plateau_region = randsample([1,2,3],1,true,[coeffs(9) coeffs(10) coeffs(11)]);
            
            switch plateau_region
                case 1
                    plateau = coeffs(6)+(coeffs(7)-coeffs(6)).*rand(1,1);
                case 2
                    plateau = coeffs(7)+(coeffs(8)-coeffs(7)).*rand(1,1);
                case 3 
                    plateau = coeffs(8)+(coeffs(9)-coeffs(8)).*rand(1,1);
                otherwise
                    plateau = coeffs(6)+(coeffs(7)-coeffs(6)).*rand(1,1);
                    warning('Unknown plateau_region, default was assigned')
            end                     
            w_d = w_n*sqrt(1-zeta^2);
            t_r = (pi-atan(w_d/(zeta*w_n)))/(w_d)
            t_inj = 50;
            [~,t_r_index] = min(abs(t-t_r));
            [~,t_inj_index] = min(abs(t-t_inj))
            if en_plateau==0
                if zeta<0.9
                    zeta = 0.9;
                end
                prof_fit = (1-(((exp(-zeta*w_n.*(t)))/sqrt(1-zeta^2)).*sin(w_d.*(t)+atan(sqrt(1-zeta^2)/zeta))))-1;
                t_np_index = find(prof_fit(1:end-1)<0 & prof_fit(2:end)>0)+1;
                prof_fit(1:t_np_index(1)) = 0;                            % wa  
                prof_fit = circshift(prof_fit,t_r_index+t_inj_index-t_np_index(1));            % wa      
                scale_factor = (coeffs(2)/max(prof_fit));
                prof_fit = scale_factor*prof_fit+basal;
                prof_fit(1:t_np_index(1)) = basal;                              
            else
                prof_fit = (1-(((exp(-zeta*w_n.*(t-t_inj-t_r)))/sqrt(1-zeta^2)).*sin(w_d.*(t-t_inj-t_r)+atan(sqrt(1-zeta^2)/zeta))));
                scale_factor = (coeffs(2)/max(prof_fit(t_inj_index+t_r_index:end)));
                prof_fit = scale_factor*prof_fit+basal;                
                prof_fit(1:t_inj_index+t_r_index) = basal;                            
            end            
            
        otherwise
            error('genProfileFit: Fit model not found')
    end