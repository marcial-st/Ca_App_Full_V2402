function [CELLSf] = featureExtraction(profiles_matrix)
    format long
    [nProf,nFrames] = size(profiles_matrix);

    CaProfsSmooth = zeros(nProf,nFrames);   %% Smooth profiles with moving average filter
    for isec = 1:1:nProf
        CaProfsSmooth(isec,:) = smooth(profiles_matrix(isec,:));
    end

    CaProfsDiff = round(diff(CaProfsSmooth')',6);  %% Ca Profs differences matrix
    maxDiff = max(max(CaProfsDiff));
    [~,tStart] = find(CaProfsDiff == maxDiff);
    tStart = tStart-4;
%     tStart = 30;

    f = waitbar(0,'Please wait...');
    
%     Structure Prealocation
    CELLSf(nProf).fStartExp = tStart;
    CELLSf(nProf).fBasalState = 999;  
    CELLSf(nProf).fMainArea = 999;
    CELLSf(nProf).fMainPeakHeight = 999;
    CELLSf(nProf).fMainPeakVal = 999;
    CELLSf(nProf).fMainPeakWidth = 999;
    CELLSf(nProf).fMainPk30 = 999;
    CELLSf(nProf).fMainPk60 = 999;
    CELLSf(nProf).fMainPk90 = 999;
    CELLSf(nProf).fTimePk30 = 999;
    CELLSf(nProf).fTimePk60 = 999;
    CELLSf(nProf).fTimePk90 = 999;
%     CELLSf(nProf).fMSError2 = 999;
%     CELLSf(nProf).fMSError3 = 999;
    CELLSf(nProf).fNumSecPeak = 999;
    CELLSf(nProf).fPeakTime = 999;
    CELLSf(nProf).fPlateau = 999;
%     CELLSf(nProf).fSlope100 = 999;
%     CELLSf(nProf).fSlope2 = 999;
%     CELLSf(nProf).fSlope3 = 999;
%     CELLSf(nProf).vSecLocs = 999;
    
    vFramesTotal = (1:1:nFrames);
    Kslope = 1000;
    Kerr = 1;
    for iP=1:1:nProf
%         CELLSf(iP).fMSError2 = 999;
%         CELLSf(iP).fMSError3 = 999;
        winBasal = tStart; % Window for Basal state 
        winPeak = 30; % Window for First Peak / Expectation: the biggest one
        thMainPeakProm = 0.1; % Threshold of Main Peak prominence [% percentage]
        thMainPeakH = 0.1;
        thSecPeakProm = 0.2;  % Threshold of Secondary Peaks prominence [% percentage]
                              % It is relative to main peak prominence: thSecPeakProm*thMainPeakProm

        CELLSf(iP).fStartExp = tStart;
                                      
        CaProfsPart2 = CaProfsSmooth(iP,tStart+winPeak:(tStart+floor((nFrames-tStart)/2)))';  % Split by 2 the remaining part of the profile 
        CaProfsPart3 = CaProfsSmooth(iP,tStart+floor((nFrames-tStart)/2):nFrames)';   % after the detected Start of experiment 
%         CaProfs100 = CaProfsSmooth(iP,nFrames-100:nFrames)';   % after the detected Start of experiment 

        [pksl,~,pkWl,pkPl] = findpeaks(CaProfsSmooth(iP,1:tStart),'MinPeakProminence',0.05,'WidthReference','halfheight'); % Get Peaks and Indices
        if~isempty(pksl) CELLSf(iP).fBnpeak = length(pksl); else CELLSf(iP).fBnpeak = 999;end
        if~isempty(pkWl) CELLSf(iP).fBwpeak = nanmean(pkWl);else CELLSf(iP).fBwpeak = 999;end
        if~isempty(pkPl) CELLSf(iP).fBppeak = nanmean(pkPl);else CELLSf(iP).fBppeak = 999;end  
        if~isempty(pkPl) CELLSf(iP).fBmxpeak = max(pkPl);else CELLSf(iP).fBmxpeak = 999;end  
        if~isempty(pkPl) CELLSf(iP).fBmnpeak = min(pkPl);else CELLSf(iP).fBmnpeak = 999;end 

        vFrames2 = vFramesTotal(tStart+winPeak:(tStart+floor((nFrames-tStart)/2)))';                % Intermediate profile section
        if isempty(vFrames2)
%             CELLSf(iP).fSlope2 = 999;
%             CELLSf(iP).fMSError2 = 999;
            CELLSf(iP).fF2c1 = 999;
            CELLSf(iP).fF2c2 = 999;
            CELLSf(iP).fF2c3 = 999;
            CELLSf(iP).fF3mean = 999;
            CELLSf(iP).fF2npeak = 999;
            CELLSf(iP).fF2wpeak = 999;
            CELLSf(iP).fF2ppeak = 999;
            CELLSf(iP).fF2mxpeak = 999;
            CELLSf(iP).fF2mnpeak = 999;  
            CELLSf(iP).fF2auc = 999;
        else
%             vF2 = [ones(length(vFrames2),1) vFrames2];
%             CELLSf(iP).fSlope2 = vF2\CaProfsPart2;
%             CELLSf(iP).fMSError2 = mean((CaProfsPart2-vF2*CELLSf(iP).fSlope2).^2);                           %% Feature Mean Square Error2
%             CELLSf(iP).fSlope2 = Kslope*CELLSf(iP).fSlope2;
%             CELLSf(iP).fMSError2 = Kerr*CELLSf(iP).fMSError2;
            t = [tStart+winPeak:1:(tStart+floor((nFrames-tStart)/2))];
            p = polyfit(t,CaProfsPart2,2);
            CELLSf(iP).fF2c1 = p(1);
            CELLSf(iP).fF2c2 = p(2);
            CELLSf(iP).fF2c3 = p(3);
            CELLSf(iP).fF2mean = mean(CaProfsPart2);
            [pksl,~,pkWl,pkPl] = findpeaks(CaProfsPart2,'MinPeakProminence',0.02,'WidthReference','halfheight'); % Get Peaks and Indices
            if~isempty(pksl) CELLSf(iP).fF2npeak = length(pksl); else CELLSf(iP).fF2npeak = 999;end
            if~isempty(pkWl) CELLSf(iP).fF2wpeak = nanmean(pkWl);else CELLSf(iP).fF2wpeak = 999;end
            if~isempty(pkPl) CELLSf(iP).fF2ppeak = nanmean(pkPl);else CELLSf(iP).fF2ppeak = 999;end
            if~isempty(pkPl) CELLSf(iP).fF2mxpeak = max(pkPl);else CELLSf(iP).fF2mxpeak = 999;end  
            if~isempty(pkPl) CELLSf(iP).fF2mnpeak = min(pkPl);else CELLSf(iP).fF2mnpeak = 999;end              
            CELLSf(iP).fF2auc = trapz(CaProfsPart2);
        end

        vFrames3 = vFramesTotal((tStart+floor((nFrames-tStart)/2)):nFrames)';               % Final profile section
        if isempty(vFrames3)
%             CELLSf(iP).fSlope3 = 999; 
%             CELLSf(iP).fMSError3 = 999;
            CELLSf(iP).fF3c1 = 999;
            CELLSf(iP).fF3c2 = 999;
            CELLSf(iP).fF3c3 = 999;
            CELLSf(iP).fF3mean = 999;
            CELLSf(iP).fF3npeak = 999;
            CELLSf(iP).fF3wpeak = 999;
            CELLSf(iP).fF3ppeak = 999;
            CELLSf(iP).fF3mxpeak = 999;
            CELLSf(iP).fF3mnpeak = 999;
            CELLSf(iP).fF3auc = 999;
        else
%             vF3 = [ones(length(vFrames3),1) vFrames3];
%             CELLSf(iP).fSlope3 = vF3\CaProfsPart3; 
%             CELLSf(iP).fMSError3 = mean((CaProfsPart3-vF3*CELLSf(iP).fSlope3).^2);              %% Feature Mean Square Error3
%             CELLSf(iP).fSlope3 = Kslope*CELLSf(iP).fSlope3;
%             CELLSf(iP).fMSError3 = Kerr*CELLSf(iP).fMSError3;
            t = [(tStart+floor((nFrames-tStart)/2)):1:nFrames];
            p = polyfit(t,CaProfsPart3,2);
            CELLSf(iP).fF3c1 = p(1);
            CELLSf(iP).fF3c2 = p(2);
            CELLSf(iP).fF3c3 = p(3);
            CELLSf(iP).fF3mean = mean(CaProfsPart3);
            [pksl,~,pkWl,pkPl] = findpeaks(CaProfsPart3,'MinPeakProminence',0.02,'WidthReference','halfheight'); % Get Peaks and Indices
            if~isempty(pksl) CELLSf(iP).fF3npeak = length(pksl); else CELLSf(iP).fF3npeak = 999;end
            if~isempty(pkWl) CELLSf(iP).fF3wpeak = nanmean(pkWl);else CELLSf(iP).fF3wpeak = 999;end
            if~isempty(pkPl) CELLSf(iP).fF3ppeak = nanmean(pkPl);else CELLSf(iP).fF3ppeak = 999;end  
            if~isempty(pkPl) CELLSf(iP).fF3mxpeak = max(pkPl);else CELLSf(iP).fF3mxpeak = 999;end  
            if~isempty(pkPl) CELLSf(iP).fF3mnpeak = min(pkPl);else CELLSf(iP).fF3mnpeak = 999;end  
            CELLSf(iP).fF3auc = trapz(CaProfsPart3);
        end
        

        vFrames100 = vFramesTotal(nFrames-100:nFrames)';                                    % 100 last samples slope
        if isempty(vFrames100)
%             CELLSf(iP).fSlope100 = 999;
            CELLSf(iP).fF100c1 = 999;
            CELLSf(iP).fF100c2 = 999;
            CELLSf(iP).fF100c3 = 999; 
            CELLSf(iP).fF100mean = 999;
        else
%             vF100 = [ones(length(vFrames100),1) vFrames100];
%             CELLSf(iP).fSlope100 = Kslope*(vF100\CaProfsSmooth(iP,nFrames-100:nFrames)');
            t = [nFrames-100:1:nFrames];
            p = polyfit(t,CaProfsSmooth(iP,nFrames-100:nFrames),2);
            CELLSf(iP).fF100c1 = p(1);
            CELLSf(iP).fF100c2 = p(2);
            CELLSf(iP).fF100c3 = p(3);
            CELLSf(iP).fF100mean = mean(CaProfsSmooth(iP,nFrames-100:nFrames));
        end        

        CELLSf(iP).fBasalState=mean(CaProfsSmooth(iP,1:winBasal));                            %% Feature Basal State
        if (mean(CaProfsPart2) > CELLSf(iP).fBasalState) 
            CELLSf(iP).fPlateau = mean(CaProfsPart2);                                         %% Feature Plateau Value
        else
            CELLSf(iP).fPlateau =0;   % TODO take zero value out of the dynamic range, 0 ->999
        end

        [pks,pkLocs,pkWidth,pkProm] = findpeaks(CaProfsSmooth(iP,:)); % Get Peaks and Indices
        [lMin,~] = islocalmin(CaProfsSmooth(iP,:));           % Get Local Minima and Prominence
%         lMinProm = nonzeros(lMinProm);
        lMinLocs = nonzeros(lMin.*vFramesTotal);    

        pkPromRelevant = find(pkProm>thMainPeakProm); % Prominence Relevant Peaks
        if isempty(pkPromRelevant)
            pkPromRelevant = find((pks-CELLSf(iP).fBasalState)>thMainPeakH);
        end
%         CELLSf(iP).fFullArea = trapz(CaProfsSmooth(iP,tStart:end));
        p_aux1 = mean(CaProfsSmooth(iP,tStart:tStart+10));
        p_aux2 = mean(CaProfsSmooth(iP,end-10:end));        
        CELLSf(iP).fSlopeFull = Kslope*(p_aux2-p_aux1)/(nFrames-tStart);

        if ~(isempty(pkPromRelevant))
            % Is there a Relevant Peak in the Window after stimulus?
            if((winBasal<=pkLocs(pkPromRelevant(1))) && (pkLocs(pkPromRelevant(1))<=(winBasal+winPeak)))
                CELLSf(iP).fMainPeakVal = CaProfsSmooth(iP,pkLocs(pkPromRelevant(1)));    %% Feature Main Peak Value
                CELLSf(iP).fMainPeakHeight = CELLSf(iP).fMainPeakVal-CELLSf(iP).fBasalState;                    %% Feature Main Peak Prominence    
                CELLSf(iP).fMainPeakWidth = pkWidth(pkPromRelevant(1));                   %% Feature Main Peak Width
                CELLSf(iP).fPeakTime = pkLocs(pkPromRelevant(1));                         %% Feature Main Peak Frame
                vSecPeak = pkProm(vFrames2(1)<pkLocs & pkLocs<vFramesTotal(end));       %% Limit indices of secondary peaks, Only 2nd section
                fNumSecPeak = (vSecPeak>=thSecPeakProm*thMainPeakProm);                 %% Filter by thresholding thSecPeakProm*thPeakProm
                vSecPeak = nonzeros(vSecPeak.*fNumSecPeak);                             %% Store secondary peaks indices
                CELLSf(iP).fNumSecPeak = sum(fNumSecPeak);                                         %% Features Number of Secondary Peaks
                if sum(fNumSecPeak)>0
                    CELLSf(iP).fSecPeakMean = mean(vSecPeak);
                else
                    CELLSf(iP).fSecPeakMean = 999;
                end                
%                 vSecLocs = zeros(size(vSecPeak));
%                 for isec=1:1:length(vSecPeak)
%                     CELLSf(iP).vSecLocs(isec) = pkLocs(pkProm == vSecPeak(isec));
%                 end

                lMinLocs = nonzeros((lMinLocs>CELLSf(iP).fPeakTime).*lMinLocs);  % Only after Main Peak
                CELLSf(iP).fMainArea = trapz(CaProfsSmooth(iP,tStart:lMinLocs(1)))-((lMinLocs(1)-tStart)*CELLSf(iP).fBasalState);
                CELLSf(iP).fMainPk90 = 0.9*CELLSf(iP).fMainPeakHeight+CELLSf(iP).fBasalState;
                CELLSf(iP).fMainPk60 = 0.6*CELLSf(iP).fMainPeakHeight+CELLSf(iP).fBasalState;
                CELLSf(iP).fMainPk30 = 0.3*CELLSf(iP).fMainPeakHeight+CELLSf(iP).fBasalState;

                auxVector=CaProfsSmooth(iP,CELLSf(iP).fPeakTime:lMinLocs(1));
                [~,fTimePk90]=min(abs(auxVector-CELLSf(iP).fMainPk90));
                [~,fTimePk60]=min(abs(auxVector-CELLSf(iP).fMainPk60));
                [~,fTimePk30]=min(abs(auxVector-CELLSf(iP).fMainPk30));

                fTimePk90 = fTimePk90+CELLSf(iP).fPeakTime;
                fTimePk60 = fTimePk60+CELLSf(iP).fPeakTime;
                fTimePk30 = fTimePk30+CELLSf(iP).fPeakTime;

                if (fTimePk90(1)>=CELLSf(iP).fPeakTime && fTimePk90(1)<=lMinLocs(1))
                    CELLSf(iP).fTimePk90 = fTimePk90(1);                
                else
                    CELLSf(iP).fTimePk90 = 999;
                    CELLSf(iP).fMainPk90 = 999;
                end

                if (fTimePk60(1)>=CELLSf(iP).fPeakTime && fTimePk60(1)<=lMinLocs(1))
                    CELLSf(iP).fTimePk60 = fTimePk60(1);
                else
                    CELLSf(iP).fTimePk60 = 999;
                    CELLSf(iP).fMainPk60 = 999;
                end

                if (fTimePk30(1)>=CELLSf(iP).fPeakTime && fTimePk30(1)<=lMinLocs(1))
                    CELLSf(iP).fTimePk30 = fTimePk30(1);
                else
                    CELLSf(iP).fTimePk30 = 999;
                    CELLSf(iP).fMainPk30 = 999;
                end                                               
            else
                CELLSf(iP).fMainPeakVal = 999;
                CELLSf(iP).fMainPeakHeight = 999;
                CELLSf(iP).fMainPeakWidth = 999;
                CELLSf(iP).fPeakTime = 999;
%                 CELLSf(iP).fNumSecPeak = 999; % NEW
                vSecPeak = pkProm(vFrames2(1)<pkLocs & pkLocs<vFramesTotal(end));       %% Limit indices of secondary peaks, Only 2nd section
                fNumSecPeak = (vSecPeak>=thSecPeakProm*thMainPeakProm);                 %% Filter by thresholding thSecPeakProm*thPeakProm
                vSecPeak = nonzeros(vSecPeak.*fNumSecPeak);                             %% Store secondary peaks indices
                CELLSf(iP).fNumSecPeak = sum(fNumSecPeak);                                         %% Features Number of Secondary Peaks
                if sum(fNumSecPeak)>0
                    CELLSf(iP).fSecPeakMean = mean(vSecPeak);
                else
                    CELLSf(iP).fSecPeakMean = 999;
                end  
%
                CELLSf(iP).fMainPk90 = 999;
                CELLSf(iP).fTimePk90 = 999;
                CELLSf(iP).fMainPk60 = 999;
                CELLSf(iP).fTimePk60 = 999;
                CELLSf(iP).fTimePk30 = 999;
                CELLSf(iP).fMainPk30 = 999;
                CELLSf(iP).fMainArea = 999;
            end
        else
                CELLSf(iP).fMainPeakVal = 999;
                CELLSf(iP).fMainPeakHeight = 999;
                CELLSf(iP).fMainPeakWidth = 999;
                CELLSf(iP).fPeakTime = 999;
%                 CELLSf(iP).fNumSecPeak = 999; % NEW
                vSecPeak = pkProm(vFrames2(1)<pkLocs & pkLocs<vFramesTotal(end));       %% Limit indices of secondary peaks, Only 2nd section
                fNumSecPeak = (vSecPeak>=thSecPeakProm*thMainPeakProm);                 %% Filter by thresholding thSecPeakProm*thPeakProm
                vSecPeak = nonzeros(vSecPeak.*fNumSecPeak);                             %% Store secondary peaks indices
                CELLSf(iP).fNumSecPeak = sum(fNumSecPeak);                                         %% Features Number of Secondary Peaks
                if sum(fNumSecPeak)>0
                    CELLSf(iP).fSecPeakMean = mean(vSecPeak);
                else
                    CELLSf(iP).fSecPeakMean = 999;
                end  
%                
                CELLSf(iP).fMainPk90 = 999;
                CELLSf(iP).fTimePk90 = 999;
                CELLSf(iP).fMainPk60 = 999;
                CELLSf(iP).fTimePk60 = 999;
                CELLSf(iP).fTimePk30 = 999;
                CELLSf(iP).fMainPk30 = 999;
                CELLSf(iP).fMainArea = 999;
        end    

        waitbar(iP/nProf,f,'Extracting Features');

    end
    close(f)     
end