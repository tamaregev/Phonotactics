
function [hf,hl] = plotmROI_phono(T,expName,effects,whichROIs,displayStats,displayIndividuals)
    %% Inputs
    % Examples:
    %   T is a data table from T=readtable('resultsfile.csv')
    %       UID                 Subject                 ROI     Effect     EffectSize
    %       ___    _________________________________    ___    ________    __________
    %
    %       18    {'018_FED_20160112d_3T2_PL2017' }     1     {'high'}      0.55596 
    %       18    {'018_FED_20160112d_3T2_PL2017' }     1     {'low' }       0.3236 
    %  
    %   effects = {'high','low'}; % should match names in Effect
    %   displayStats = true;% or false
    %   displayIndividuals = true;%or false
    %   whichROIs = [2:6];%choose manually

    
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
    %% definitions
    plotT = T(ismember(T.ROI,whichROIs) & ismember(T.Effect,effects),:); 
    nSubj=numel(unique(plotT.Subject));
    EffectSize = reshape(plotT.EffectSize,[numel(effects),nSubj,numel(whichROIs)]);
    meanE = squeeze(mean(EffectSize,2));
    stderr = squeeze(std(EffectSize,0,2)./sqrt(size(EffectSize,2)));
    switch expName
        case 'langlocSN'
            colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
        case 'SWJNaud'
            colors = [0 0.1470 0.7410; 0.1235 0.4745 0.8667; 0.6196 0.6902 0.9333; 0.8500 0.3250 0.0980];
        case 'Biling'
            colors = [1 0.3 0.1; 0.8500 0.3250 0.0980];
        case 'ParamNew'
            order = 10:-1:1;%reorder conditions from S-S1-4, W, W1-4
            EffectSize = EffectSize(order,:,:);
            meanE = squeeze(mean(EffectSize,2));
            stderr = squeeze(std(EffectSize,0,2)./sqrt(size(EffectSize,2)));
            colors = [187 33 6; 204 77 57; 221 121 107; 237 164 158; 254 208 208; 6 12 186; 57 67 203; 108 121 221; 158 176 238; 209 230 255]./255;
       case 'Nhood'
            colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];        
    end
    %% plot
    addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/downloaded')%this should include barwitherr function downloaded from matheorks
    
    %calc significance
%     p = nan(length(whichROIs),1); H = nan(length(whichROIs),1); CI = cell(length(whichROIs),1); stats = cell(length(whichROIs),1);
%     for iroi = 1:length(whichROIs)
%         ROI = whichROIs(iroi);
%         [p(iroi),H(iroi),CI{iroi},stats{iroi}] = significance_mROI(T_ind,contrast1,contrast2, ROI);
%     end
%     mROIsig=table(whichROIs,p,H,CI,stats);
%     

   
    hf=figure;set(hf,'Position',[100 100 1300 700])
    [h, hE]=barwitherr(stderr',1:length(whichROIs),meanE'); 
    for ib=1:length(h)
        h(ib).FaceColor=colors(ib,:);
        h(ib).LineWidth=2;
    end
    set(hE,'Linewidth',2)

    xlabel('ROI')
    ylabel('Effect size')

    set(gca,'xtick',1:length(whichROIs))
    set(gca,'xticklabels',whichROIs)
    
    hold on
    if displayIndividuals
        for ib=1:numel(h)
            XData = h(ib).XData+h(ib).XOffset;
            YData = h(ib).YData;
            for iroi=1:numel(XData)
                switch expName
                    case 'Nhood'
                        ind = squeeze(EffectSize(ib,:,iroi));
                    case 'ParamNew'
                        ind = squeeze(EffectSize(ib,:,iroi));
                    case 'langlocSN'
                        ind = squeeze(EffectSize(ib,:,iroi));
                    case 'SWJNaud'
                        ind = squeeze(EffectSize(ib,:,iroi));
                     case 'Biling'
                        ind = squeeze(EffectSize(ib,:,iroi));
                   
                end
                xx=repmat(XData(iroi),size(ind));
   %             hs=plot((xx+(rand(size(xx))-0.5)/10),ind,'o','MarkerFaceColor',colors(ib,:),'MarkerEdgeColor',colors(ib,:).*0.5,'MarkerSize',6,'Linewidth',0.1,'HandleVisibility','off');
                switch expName
                    case 'langlocSN'
                        hs=scatter((xx+(rand(size(xx))-0.5)/10),ind,30,colors(ib,:).*0.8,'filled','HandleVisibility','off');
                        hs.MarkerFaceAlpha=0.1;
                    case 'SWJNaud'
                        hs=scatter((xx+(rand(size(xx))-0.5)/100),ind,25,colors(ib,:).*0.8,'filled','HandleVisibility','off');
                        hs.MarkerFaceAlpha=0.2;
                    case 'ParamNew'    
%                        hs=scatter((xx+(rand(size(xx))-0.5)/100),ind,25,colors(ib,:).*0.8,'filled','MarkerEdgeColor',colors(ib,:).*0.5,'HandleVisibility','off');
                        hs=scatter((xx+(rand(size(xx))-0.5)/100),ind,25,colors(ib,:).*0.8,'filled','HandleVisibility','off');
                        hs.MarkerFaceAlpha=0.2;
                    case 'Nhood'
                        hs=scatter((xx+(rand(size(xx))-0.5)/10),ind,30,colors(ib,:).*0.8,'filled','HandleVisibility','off');
                        hs.MarkerFaceAlpha=0.2;
                    case 'Biling'
                        hs=scatter((xx+(rand(size(xx))-0.5)/10),ind,30,colors(ib,:).*0.8,'filled','HandleVisibility','off');
                        hs.MarkerFaceAlpha=0.2;
               
                end
                hE2=errorbar(XData(iroi),YData(iroi),stderr(ib,iroi),'Color','k');
                set(hE2,'linewidth',2)
            end
        end
        ylim([min(min(min(EffectSize))), max(max(max(EffectSize)))])
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
    
    hl=legend(strrep(effects,'_',' '),'Location','neo');
    set(hl,'fontsize',20)
    set(gca,'fontsize',16)
   
end

