
function [h,hE,hl] = plotmROI_phono_allexp(T,effects,whichROIs,ROIstring,whicherr,avgfROIs,displayStats,displayIndividuals)
    %% Inputs
    % Examples:
    %   T is a data table from T=readtable('resultsfile.csv')
    %      UID                Subject                 ROI       Effect        EffectSize    ExpNum
    %      ___    ________________________________    ___    _____________    __________    ______
    %
    %      18    {'018_FED_20160112d_3T2_PL2017' }     1       {'high'}        0.55596       1 
    %      18    {'018_FED_20160112d_3T2_PL2017' }     1       {'low' }        0.3236        1     
    %  
    %   effects = {{'N langloc'},{'W','W1','W2','W3','W4'},{'high','low'}}; % should match names in Effect,
    %              ordered due to ExpNum
    %   whichROIs = 1:6;
    %   ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
    %   whicherr = 'stderr' or 'confidence';
    %   avgfROIs = true or false. If true - plots the avg of all fROIs
    % . displayStats = true;% or false (TODO)
    %   displayIndividuals = true;%or false
    % 
    %% manage inputs
    if ~exist('effects','var')
        effects = unique(T.Effect);%display as text over figure
    end
    
    if ~exist('displayStats','var')
        displayStats = false;%display as text over figure
    end
    if ~exist('displayIndividuals','var')
        displayIndividuals = true;%display as text over figure
    end

    if ~exist('whichROIs','var')
        whichROIs = 1:6;%select after reading data due to several criteria
    end
    if ~exist('avgfROIs','var')
        avgfROIs = false;%select after reading data due to several criteria
    end
    if length(whichROIs)==1
        if avgfROIs
            error('Inputs not compatible: averaging fROIs is not compatible with passing only 1 fROI.')
        end
    end
    %% calc
     addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')
     
    clear allEffects
    allEffects = cell(0,0);
    for ii=1:length(effects)
        allEffects = [allEffects,effects{ii}];
    end
    whichExpPerEffect = nan(size(allEffects));
   
    colors = [255 255 255; 150 150 150; 100 100 100; 6 12 186; 57 67 203; 108 121 221; 158 176 238; 209 230 255; 247 37 133;251 157 199]./255;
    
    nexp = numel(effects);

    EffectSize = cell(nexp,1);
    meanE_all = nan(numel(allEffects),numel(whichROIs));
    err_all = nan(numel(allEffects),numel(whichROIs));
    meanE_allfROIs = nan(numel(allEffects),1);
    err_allfROIs = nan(numel(allEffects),1);
    
    ieffect = 1;
    for iexp=1:nexp
        plotT = T(ismember(T.ROI,whichROIs) & ismember(T.Effect,allEffects) & T.ExpNum==iexp ,:); 
        allSubj = unique(plotT.UID);
        nSubj=numel(allSubj);
        EffectSize{iexp}=nan(numel(effects{iexp}),nSubj,numel(whichROIs));
        for ii=1:height(plotT)
            ief=find(ismember(effects{iexp},plotT.Effect(ii)));
            isub=find(ismember(allSubj,plotT.UID(ii)));
            iroi=find(ismember(whichROIs,plotT.ROI(ii)));
            EffectSize{iexp}(ief,isub,iroi)=plotT.EffectSize(ii);    
        end
        meanE{iexp} = mean(EffectSize{iexp},2);
        
        stderr{iexp} = (std(EffectSize{iexp},0,2)./sqrt(size(EffectSize{iexp},2)));
        conf{iexp} = nan(size(stderr{iexp}));
        for iroi=1:length(whichROIs)
            es=EffectSize{iexp}(:,:,iroi)';
            conf{iexp}(:,:,iroi) = Confidence(es);
        end
        meanE_all(ieffect:ieffect + numel(effects{iexp}) - 1 ,:) = meanE{iexp};
        switch whicherr
            case 'confidence'
                err_all(ieffect:ieffect + numel(effects{iexp}) - 1 ,:) = conf{iexp};
            case 'stderr'
                err_all(ieffect:ieffect + numel(effects{iexp}) - 1 ,:) = stderr{iexp};
        end
        
        if length(whichROIs)>1 %this is calculated for the avgfROIs option
            meanE_perfROI{iexp} = mean(EffectSize{iexp},3);
            switch whicherr
                case 'confidence'
                    errs = Confidence(meanE_perfROI{iexp}');
                case 'stderr'
                    errs = std(meanE_perfROI{iexp},0,2)./sqrt(size(meanE_perfROI{iexp},2)); 
            end
            err_allfROIs(ieffect:ieffect + numel(effects{iexp}) - 1 ) = errs;
        end
        
        whichExpPerEffect(ieffect:ieffect + numel(effects{iexp}) - 1) = ones(size(effects{iexp}))*iexp;
        ieffect = ieffect + numel(effects{iexp});
    end    
    meanE_allfROIs = mean(meanE_all,2);
    
    %% plot
    if avgfROIs
        for iexp=1:nexp
            EffectSize{iexp} = mean(EffectSize{iexp},3);
        end   
        Plot(meanE_allfROIs,err_allfROIs,EffectSize)
    else
        Plot(meanE_all,err_all,EffectSize)
    end
    %% nested functions
    function Plot(means,errs,subjData)
    
        addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/downloaded')%this should include barwitherr function downloaded from matheorks

        %calc significance
    %     p = nan(length(whichROIs),1); H = nan(length(whichROIs),1); CI = cell(length(whichROIs),1); stats = cell(length(whichROIs),1);
    %     for iroi = 1:length(whichROIs)
    %         ROI = whichROIs(iroi);
    %         [p(iroi),H(iroi),CI{iroi},stats{iroi}] = significance_mROI(T_ind,contrast1,contrast2, ROI);
    %     end
    %     mROIsig=table(whichROIs,p,H,CI,stats);
   
        %insert nans to create spaces between the 3 experiments:
        nSpaces = nexp-1;%add due to number of experiments
        meanE_all_wnan = nan(size(means,1)+nSpaces,size(means,2));
        err_all_wnan = nan(size(means,1)+nSpaces,size(means,2));
        for ie=1:nexp
            indfrom_wn = 1;
            indto_wn = numel(effects{1});
            indfrom = 1;
            indto = numel(effects{1});
            if ie>1
                for i=2:ie
                    indfrom_wn = indfrom_wn + numel(effects{i-1}) + 1;
                    indto_wn = indto_wn + numel(effects{i}) + 1;
                    
                    indfrom = indfrom + numel(effects{i-1});
                    indto = indto + numel(effects{i});

                end
            end
            meanE_all_wnan(indfrom_wn:indto_wn,:) = means(indfrom:indto,:);
            err_all_wnan(indfrom_wn:indto_wn,:) = errs(indfrom:indto,:);
        end

        ibMap=find(~isnan(meanE_all_wnan(:,1)));

        if size(means,2)>1 %plot several fROIs
            [h, hE]=barwitherr(err_all_wnan',1:length(whichROIs),meanE_all_wnan'); 
            for ib=1:length(h)
                if find(ibMap==ib)
                    h(ib).FaceColor=colors(find(ibMap==ib),:);
                    h(ib).LineWidth=1;
                else
                    h(ib).HandleVisibility = 'off';
                end
            end
            set(hE,'Linewidth',1)
            set(gca,'xtick',1:length(whichROIs))
            xticklabels(gca,ROIstring(whichROIs))
            xlabel('ROI')
            
        else %plot 1 or avg fROI
            for ib=1:length(err_all_wnan)
                if find(ibMap==ib)
                    h(ib) = bar(ib,meanE_all_wnan(ib)); 
                    h(ib).FaceColor=colors(find(ibMap==ib),:);
                    hold on
                    hE(ib) = errorbar(ib,meanE_all_wnan(ib),err_all_wnan(ib));
                    set(hE(ib),'Color','k','HandleVisibility','off')
                    if avgfROIs
                        h(ib).LineWidth=1;
                        set(hE(ib),'linewidth',1)
                    else
                        h(ib).LineWidth=1;
                        set(hE(ib),'linewidth',1)
                    end
                   
                end
            end
            set(gca,'xtick',[1 3 5 9 13.5])
            ExpString = {'Exp 1a','Exp 1b','Exp 1c','Exp 2','Exp 3'};
            Nstring = {'N=605','N=12','N=13','N=16','N=14'};
            labelArray = [ExpString;Nstring];
            tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
            set(gca,'xticklabels',tickLabels)
            if avgfROIs
                title('Averaged across fROIs')
            else
                title(ROIstring(whichROIs))
            end
        end
        set(gca,'fontsize',18)
        ylabel('% BOLD signal change','fontsize',16)
        

        hold on
        if displayIndividuals
            if size(means,2)>1 %plot several fROIs
                circleSize = [20 30];
                spreadScale = [30 30];
                alphaLevel = [0.02 0.2];
            elseif avgfROIs
                circleSize = [50 70];
                spreadScale = [3 5];    
                alphaLevel = [0.02 0.2];
            else %plot 1 fROI
                circleSize = [7 14];
                spreadScale = [2 3];
                alphaLevel = [0.02 0.2];
            end
            for ib=1:numel(h)
                if find(ibMap==ib)
                    iib = find(ibMap==ib);
                    XData = h(ib).XData+h(ib).XOffset;
                    YData = h(ib).YData;
                    iexp = whichExpPerEffect(iib);
                    ibPerExp = iib;
                    if iexp>1
                        for i=2:iexp
                            ibPerExp = ibPerExp - numel(effects{i-1});
                        end
                    end
%                     
%                     switch iexp
%                         case 1
%                             ibPerExp = iib;
%                         case 2 
%                             ibPerExp = iib - numel(effects{1});
%                         case 3
%                             ibPerExp = iib - numel(effects{1}) - numel(effects{2});
%                     end
                    
                    for iroi=1:numel(XData)

                        ind = squeeze(subjData{iexp}(ibPerExp,:,iroi));

                        xx=repmat(XData(iroi),size(ind));
           %             hs=plot((xx+(rand(size(xx))-0.5)/10),ind,'o','MarkerFaceColor',colors(ib,:),'MarkerEdgeColor',colors(ib,:).*0.5,'MarkerSize',6,'Linewidth',0.1,'HandleVisibility','off');
                        switch iexp
                            case 1
                                hs=scatter((xx+(rand(size(xx))-0.5)/spreadScale(1)),ind,circleSize(1),'k','filled','HandleVisibility','off');
                                hs.MarkerFaceAlpha=alphaLevel(1);
                            otherwise  
        %                        hs=scatter((xx+(rand(size(xx))-0.5)/100),ind,25,colors(ib,:).*0.8,'filled','MarkerEdgeColor',colors(ib,:).*0.5,'HandleVisibility','off');
                                hs=scatter((xx+(rand(size(xx))-0.5)/spreadScale(2)),ind,circleSize(2),colors(iib,:).*0.8,'filled','HandleVisibility','off');
                                hs.MarkerFaceAlpha=alphaLevel(2);
                                %hs=scatter((xx+(rand(size(xx))-0.5)/spreadScale(2)),ind,circleSize(2),colors(iib,:).*0.8,'filled','HandleVisibility','off');
                                %hs.MarkerFaceAlpha=alphaLevel(2);
                        end
    %                hE2=errorbar(XData(iroi),YData(iroi),stderr(ib,iroi),'Color','k');
    %                set(hE2,'linewidth',2)
                    end
                end
            end
            allmins = nan(size(effects));
            for ie = 1:numel(effects)
                allmins(ie) = min(min(min(subjData{ie})));
                allmaxs(ie) = max(max(max(subjData{ie})));
            end
            %minind = min([min(min(min(subjData{1}))) min(min(min(subjData{2}))) min(min(min(subjData{3})))]);
            minind = min(allmins);
%            maxind = max([max(max(max(subjData{1}))) max(max(max(subjData{2}))) max(max(max(subjData{3})))]);
            maxind = max(allmaxs);

            ylim([minind, maxind])
            %ylim([-2, 4])
        end


        if displayStats
            if displayIndividuals
                dy=-2;dx=0;fs=10;
            else
                dy=-1;dx=0;fs=10;
            end

            for iroi=1:length(whichROIs)
                roi=whichROIs(iroi);
                %text(iroi,max(mean(iroi,:))+10*dy,num2str(roi))
                str = {['p mROI=' num2str(mROIsig.p(iroi),2)],[ 'overlap: ' num2str(T_GSS.inter_subjectOverlap(roi),2)],[ 'p GSS = ' num2str(T_GSS.p(roi),2)],['p fdr = ' num2str(T_GSS.p_fdr(roi),2)]};
                if mROIsig.H(iroi)
                    c = [1 0 0];
                elseif mROIsig.p(iroi) <= 0.1
                    c = [0 0 1];
                else
                    c = [0 0 0];
                end
                if numel(whichROIs) > 1
                    t=text(iroi-dx,dy,str,'fontsize',fs,'fontweight','bold','horizontalalignment','center','Color',c);
                else
                    t=text(numel(whichEffects)/2+0.5,dy,str,'fontsize',fs,'fontweight','bold','horizontalalignment','center','Color',c);
                end       
             end
             if ~displayIndividuals
                cy=get(gca,'ylim');
                ylim([2*dy, max(max(mean+stderr))])
             else
                cy=get(gca,'ylim');
                if cy(1) > 2*dy
                    ylim([2*dy, cy(2)])
                end
             end

        end

        hl=legend(strrep(allEffects,'_',' '),'Location','neo');

    %     for iexp = 1:length(effects)
    %         ieffects=find(whichExpPerEffect==iexp);
    %         hl(iexp)=legend([h(ibMap(ieffects))],strrep(allEffects(ieffects),'_',' '),'Location','neo');
    %         pos=get(hl(iexp),'Position');
    %         pos(2) = pos(2) - (iexp-1)/10;
    %         hold on
    %         set(hl(iexp),'Position',pos);
    %     end    
        set(hl,'fontsize',20)
    end
end

