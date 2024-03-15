function [pks_all,pks_loc_all,aux_all_p,aux_all_n,filt_a_p,filt_a_n]=th_fine_approach_v1(profiles_smooth,th_p)
    [np,nf] = size(profiles_smooth);
    REPORT=false;
    deltaMx=3;
    delta=7;%GT to estimated max
    lfMostMx=5;
%     PromTh= 0.01825;%0.0175  0.01250
    timesProm=4;
    maxWidth= 24;
    minSep=7;
%     nCells= length(profiles_matrix(:,1));
    samples=(1:length(profiles_smooth(1,:)));
    nSamples=length(samples);
    OddProf=0;
    for j=1:1:np%[25 45 60 75 81 91 103]%1:nCells %[60 81 99]%[105 109 119 121 123 127 128 160 163 180]%1:nCells %106 110 [7 9 40 42 89 91]%[2 30 48 51 54 61 111]%
        %==============Analyze smooth profile by Differences=============
        mDif=diff(profiles_smooth(j,:));
        vDiff=double(mDif>0.0).*mDif+1.0;
        didx=1;
        nMin=0;
        nMax=0;
        descIdx=0;
        minCONT=7;
        nSamples=length(samples);
        
        
        %         --------------------------
        lDiff= double(vDiff > 1)*1.25 + double(vDiff ==1);
        %----------------------
        lastVal=0;
        cPos=0;
        cNeg=0;
        PosIntval=[];
        NegIntval=[];
        for iFr=1: length(lDiff)
            if (lastVal == lDiff(iFr))
                continue;
            else
                if (lDiff(iFr)==1.25)
                    cPos=cPos+1;
                    PosIntval(cPos,1)=iFr;
                    if (cNeg >0)
                        NegIntval(cNeg,2)=iFr-1;
                    end
                else
                    cNeg=cNeg+1;
                    NegIntval(cNeg,1)=iFr;
                    if (cPos >0)
                        PosIntval(cPos,2)=iFr-1;
                    end
                end
                lastVal=lDiff(iFr);
            end
        end
        posDiff=PosIntval(:,2)-PosIntval(:,1)+1;
        negDiff=NegIntval(:,2)-NegIntval(:,1)+1;
        nDiff=lDiff;
        for didx=1: length(negDiff)
            if (negDiff(didx)== 1)
                nDiff(didx)= 1.25;
            end
            %             if (posDiff(didx)== 1)
            %                nDiff(didx)= 1.0;
            %             end
        end
        %----------------------
        lastVal=0;
        cPos=0;
        cNeg=0;
        sPosIntval=[];
        sNegIntval=[];
        for iFr=1: length(nDiff)
            if (lastVal == nDiff(iFr))
                continue;
            else
                if (nDiff(iFr)==1.25)
                    cPos=cPos+1;
                    sPosIntval(cPos,1)=iFr;
                    if (cNeg >0)
                        sNegIntval(cNeg,2)=iFr-1;
                    end
                else
                    cNeg=cNeg+1;
                    sNegIntval(cNeg,1)=iFr;
                    if (cPos >0)
                        sPosIntval(cPos,2)=iFr-1;
                    end
                end
                lastVal=nDiff(iFr);
            end
        end
        snegDiff=sNegIntval(:,2)-sNegIntval(:,1)+1;
        for didx=1: length(snegDiff)
            if (snegDiff(didx)== 1)
                nDiff(sNegIntval(didx,2))= 1.25;
            end
        end
        %+++++++++++++++++++++++++++++
        lastVal=0;
        cPos=0;
        cNeg=0;
        sPosIntval=[];
        sNegIntval=[];
        for iFr=1: length(nDiff)
            if (lastVal == nDiff(iFr))
                continue;
            else
                if (nDiff(iFr)==1.25)
                    cPos=cPos+1;
                    sPosIntval(cPos,1)=iFr;
                    if (cNeg >0)
                        sNegIntval(cNeg,2)=iFr-1;
                    end
                else
                    cNeg=cNeg+1;
                    sNegIntval(cNeg,1)=iFr;
                    if (cPos >0)
                        sPosIntval(cPos,2)=iFr-1;
                    end
                end
                lastVal=nDiff(iFr);
            end
        end
        
        %________________
        sposDiff=sPosIntval(:,2)-sPosIntval(:,1)+1;
        for didx=1: length(sposDiff)
            if (sposDiff(didx)<= 2)
                if (sPosIntval(didx,2) >0)
                    nDiff(sPosIntval(didx,2))= 1.0;
                else
                    nDiff(sPosIntval(didx,1))= 1.0;
                end
                if (sposDiff(didx)== 2)
                    nDiff(sPosIntval(didx,2)-1)= 1.0;
                end
            end
        end
        %+++++++++++++++++++++++++++++
        lastVal=0;
        cPos=0;
        cNeg=0;
        sPosIntval=[];
        sNegIntval=[];
        for iFr=1: length(nDiff)
            if (lastVal == nDiff(iFr))
                continue;
            else
                if (nDiff(iFr)==1.25)
                    cPos=cPos+1;
                    sPosIntval(cPos,1)=iFr;
                    if (cNeg >0)
                        sNegIntval(cNeg,2)=iFr-1;
                    end
                else
                    cNeg=cNeg+1;
                    sNegIntval(cNeg,1)=iFr;
                    if (cPos >0)
                        sPosIntval(cPos,2)=iFr-1;
                    end
                end
                lastVal=nDiff(iFr);
            end
        end
        
        
        sposDiff=sPosIntval(:,2)-sPosIntval(:,1)+1;
        %-------------------------------
        %--------Obtaining Maxs and then mins between Maxs------------------
        %--- First locate nearest max to estimated Max with +- deltaMx---
        %--- Using sposDiff right limit------
        mMaxs=[];
        % mMaxs(1)=sPosIntval(1,2); % checar
        for iMx=1: length(sPosIntval(:,1))
            if (sPosIntval(iMx,2) > deltaMx && sPosIntval(iMx,2)+ deltaMx < nSamples)
                
                [null, mMaxs(iMx)]= max(profiles_smooth(j,sPosIntval(iMx,2)-deltaMx:sPosIntval(iMx,2)+deltaMx));
                mMaxs(iMx)= mMaxs(iMx)+sPosIntval(iMx,2)-deltaMx-1;
            else
                mMaxs(1)=sPosIntval(1,2);
            end
        end
        
        % -- Get mins between maxs---------
        kMins=[];
        nonMaxs=[];
        iMin=0;
        if (mMaxs(1) > lfMostMx)
            [null,kMins]=min(profiles_smooth(j,1:mMaxs(1)-1));
        end
        for iMx=1:length(mMaxs)-1
            [null,  LMin]=min(profiles_smooth(j,mMaxs(iMx)+1:mMaxs(iMx+1)-1));
            LMin=LMin + mMaxs(iMx);
            %             rExp{i}.Activities(j,iMx)= {[mMins(LMin), mMaxs(iMx), mMins(RMin)]};
            if (profiles_smooth(j,mMaxs(iMx+1))-profiles_smooth(j,LMin) >= th_p)
                kMins = [kMins LMin ];
            else
                nonMaxs=[nonMaxs  (iMx+1)];
            end
        end
        if (~isempty(kMins)  && kMins(1) < mMaxs(1) && profiles_smooth(j,mMaxs(1))-profiles_smooth(j,1) < th_p)
            nonMaxs= [1 nonMaxs];
            kMins(1)=[];
        end
        aux_all_n = nonMaxs;
        mMaxs(nonMaxs)=[];
        
%                 gtPks=[];
%         if (~isempty(spks(j).pks))
%         gtPks=round(spks(j).pks(:,1)/3);
%         kMatch=0;
%         
%         for ii=1: length(mMaxs)
%             gtMatch=find( gtPks >= mMaxs(ii)-delta & gtPks <= mMaxs(ii) +delta,1,'first');
%             if (~isempty(gtMatch))
%                 kMatch= kMatch +1;
%             end
%         end
%         end
% 
%         if (~isempty(gtPks))
%             pctMatch= 100*kMatch/length(gtPks);
%         else
%             pctMatch=100;
%         end       
        pks_all = [profiles_smooth(mMaxs) profiles_smooth(nonMaxs)];
        pks_loc_all = [mMaxs nonMaxs];
        aux_all_p = mMaxs;
        filt_a_p = false(1,length(pks_all));
        filt_a_n = false(1,length(pks_all));
        filt_a_p(1:length(mMaxs)) = true;
        filt_a_n(length(mMaxs)+1:end) = true;
        
%         %............................
%         %++++++++++++++
%         extras=find(gtPks > nSamples);
%         if (~isempty(extras))
%             gtPks(extras)=nSamples;
%         end
%         figure;
%         if (~isempty(gtPks) && ~isempty(mMaxs))
%             if (~REPORT)
%                 p=plot(samples*3, profiles_smooth(j,:),...
%                     mMaxs*3,profiles_smooth(j,mMaxs),'r<',...
%                     kMins*3,profiles_smooth(j,kMins),'ko',...
%                     gtPks*3, profiles_smooth(j,gtPks),'gp',...
%                     samples(1:end-1)*3, vDiff,...
%                     samples(1:end-1)*3, lDiff,...
%                     samples(1:end-1)*3, nDiff);
%                 yticks(0.75:0.05:1.55);
%                 ylim([0.75 1.55]);
%                 p(1).LineWidth=2;
%                 p(1).Color='b';
%                 p(6).Color='r';
%                 if (~isempty(kMins))
%                     p(3).LineWidth=2;
%                 end
%                 if (~isempty(gtPks))
%                     p(4).LineWidth=2;
%                     p(7).LineWidth=2;
%                 end
%             else
%                 p=plot(samples*3, profiles_smooth(j,:),...
%                     mMaxs*3,profiles_smooth(j,mMaxs),'r<',...
%                     kMins*3,profiles_smooth(j,kMins),'ko',...
%                     gtPks*3, profiles_smooth(j,gtPks),'gp');
%                 yticks(0.75:0.05:1.55);
%                 ylim([0.75 1.55]);
%                 p(1).LineWidth=2;
%                 p(1).Color='b';
%                 if (~isempty(kMins))
%                     p(3).LineWidth=2;
%                 end
%                 if (~isempty(gtPks))
%                     p(4).LineWidth=2;
%                 end
%             end
%         else
%             p=plot(samples*3, profiles_smooth(j,:));
%             yticks(0.75:0.05:1.55);
%             ylim([0.75 1.55]);
%             p(1).LineWidth=2;
%         end
%         title([ Exp{i} '   Profile# ' int2str(j)  ' Prom=' num2str(PromTh) ' %Atos=' num2str(pctMatch)]);
%         %------- 2 images per page--------
%         if (REPORT)
%             fig = Figure();
%             if (OddProf ==0)
%                 imTop0 = Image(getSnapshotImage(fig,rpt));
%                 imTop0.Height='4.5in';
%                 delete(gcf);
%                 OddProf= OddProf+1;
%             else
%                 imTop1 = Image(getSnapshotImage(fig,rpt));
%                 imTop1.Height='4.5in';
%                 delete(gcf);
%                 OddProf= OddProf+1;
%             end
%             if (OddProf==2)
%                 t = Table({imTop0;imTop1});
%                 append(rpt,t);
%                 OddProf=0;
%             end
%             if (j == nCells && OddProf==1)
%                 append(rpt,imTop0);
%                 OddProf=0;
%             end
%         end
    end