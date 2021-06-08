%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir = [root_dir 'figures/'];
expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'};
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')
%% read and plot each exp separately
for iexp=1:length(expNames)
    expName = expNames{iexp};
    T=readtable([results_dir expName '.csv']);
    switch expNames{iexp}
        case 'langlocSN'    
            ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
            [h,hl]=plotmROI_phono(T,expName,{'S','N'},1:12);
        case 'SWJNaud'
            effects={'S','W','J','N'};        
            ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
            [h,hl]=plotmROI_phono(T,expName,effects,1:12);
        case 'Biling'
            effects = {'ENGNW','CHNW'};        
            ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
            [h,hl]=plotmROI_phono(T,expName,effects,1:12);
        
        case 'ParamNew'
            effects={'S','S1','S2','S3','S4','W','W1','W2','W3','W4'};        
            ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
            [h,hl]=plotmROI_phono(T,expName,effects,1:12);
        case 'Nhood'
            ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
            [h,hl]=plotmROI_phono(T,expName,{'high','low'},1:12);
    end
    title(expName,'fontsize',24);
    xticklabels(gca,ROIstring)
    set(gca,'fontsize',16)
    grid on
    set(h,'Position',[100 100 2000 500])
    saveas(gcf,[figures_dir expName],'fig')
    saveas(gcf,[figures_dir expName],'png')
end
%% individuals Nhood

expName = 'Nhood';
T=readtable([results_dir expName '.csv']);
UIDs=unique(T.UID);
for uid=UIDs'
   h=figure;set(h,'Position',[100 100 2000 500])
   bar(reshape(T.EffectSize(T.UID==uid,:),2,12)')   
   set(gca,'xticklabels',ROIstring)
   set(gca,'fontsize',12)
   title(['UID = ' num2str(uid) ', Nhood'])
   legend({'high','low'})
    saveas(gcf,[figures_dir expName '_UID' num2str(uid)],'png')
end

%% Make and save table of all 5 Exp together
clear T
for iexp=1:length(expNames)
    expName = expNames{iexp};
    expN = expNs{iexp};
    T{iexp}=readtable([results_dir expN filesep expName '.csv']);
    ExpNum=iexp*ones(height(T{iexp}),1);
    ExpName = cell(size(ExpNum));
    ExpName(:,:) = {expName}; 
    T{iexp}=[T{iexp} table(ExpNum) table(ExpName)];
    switch expName
       case 'langlocSN'    
           T{iexp}.Effect(strcmp(T{iexp}.Effect,{'S'}))={'S langloc'};
           T{iexp}.Effect(strcmp(T{iexp}.Effect,{'N'}))={'N langloc'};
       case 'SWJNaud'
           T{iexp}.Effect(strcmp(T{iexp}.Effect,{'S'}))={'S aud'};
           T{iexp}.Effect(strcmp(T{iexp}.Effect,{'N'}))={'N aud'};
           T{iexp}.Effect(strcmp(T{iexp}.Effect,{'W'}))={'W aud'};
           T{iexp}.Effect(strcmp(T{iexp}.Effect,{'J'}))={'J aud'};
       
           T{iexp}=T{iexp}(:,[2,1,3,5,6,7]);
           UID=nan(height(T{iexp}),1);
           for ii=1:height(T{iexp})
               UID(ii)=str2num(T{iexp}.Subject{ii}(1:3));
           end
           T{iexp}=[table(UID), T{iexp}];
        case 'Biling'
       
           T{iexp}=T{iexp}(:,[2,1,3,5,6,7]);
           UID=nan(height(T{iexp}),1);
           for ii=1:height(T{iexp})
               UID(ii)=str2num(T{iexp}.Subject{ii}(1:3));
           end
          T{iexp}=[table(UID), T{iexp}];
       case 'ParamNew'
           T{iexp}=T{iexp}(:,[2,1,3,5,6,7]);
           UID=nan(height(T{iexp}),1);
           for ii=1:height(T{iexp})
               UID(ii)=str2num(T{iexp}.Subject{ii}(1:3));
           end
           T{iexp}=[table(UID), T{iexp}];
       case 'Nhood'
           T{iexp}=T{iexp}(:,[1:4,6:8]);
    end
    ExpN = cell(size(ExpNum));
    ExpN(:,:) = {expN}; 
    T{iexp}=[T{iexp} table(ExpN)];
end
Tb = [T{1}; T{2}; T{3}; T{4}; T{5}];

save([results_dir 'Tb'],'Tb')
%% plot all Exp together
tag = 'LH';
hf=figure;set(hf,'Position',[100 100 1500 500])

switch tag
    case 'LH'
        whichROIs=1:6;
    case 'RH'
        whichROIs=7:12;
    case 'all'
        whichROIs=1:12;
end
       
effects = {{'N langloc'},{'N aud'},{'ENGNW'},{'W','W1','W2','W3','W4'},{'high','low'}};
ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
[hf,hl] = plotmROI_phono_allexp(Tb,effects,whichROIs,ROIstring,'confidence');
title(tag)
ylim([-2 6])
saveas(gcf,[figures_dir 'allExp_' tag],'fig')
saveas(gcf,[figures_dir 'allExp_' tag],'png')
%% .. 1 fROIs
whichROIs=2;
hf=figure;set(hf,'Position',[100 100 800 500])
ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
[h,hE,hL] = plotmROI_phono_allexp(Tb,effects,whichROIs,ROIstring,'confidence');
%% avg fROI
whichROIs=1:6;
avgfROIs = true;
hf=figure;set(hf,'Position',[100 100 800 500])
[h,hE,hL] = plotmROI_phono_allexp(Tb,effects,whichROIs,ROIstring,'confidence',avgfROIs);
%% plot combination
load([results_dir 'Tb'])

tag='comb';
effects = {{'N langloc'},{'N aud'},{'ENGNW'},{'W','W1','W2','W3','W4'},{'high','low'}};
ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
rois = 7:11;

hf=figure;set(hf,'Position',[100 100 800 800])
subplot(4,6,[1:6])
whichROIs=rois;
avgfROIs = true;
[h,hE,hL] = plotmROI_phono_allexp(Tb,effects,whichROIs,ROIstring,'stderr',avgfROIs);
set(hL,'fontsize',16)
set(gca,'Position',[0.15 0.65 0.65 0.3])
set(gca,'Ytick',[0,1,2])
ylabel('% BOLD signal change','fontsize',16)
        
avgfROIs = false;
spi = {13:14,15:16,17:18,20:21,22:23};
for ir=1:length(rois)
    whichROI=rois(ir);
    subplot(4,6,spi{ir})
    %[h,hE,hL] = plotmROI_phono_allexp(Tb,effects,whichROIs,ROIstring,'confidence');
    [h,hE,hL] = plotmROI_phono_allexp(Tb,effects,whichROI,ROIstring,'stderr');
    legend off
    ylim([-2 3])
    if whichROI~=6
        set(gca,'XTickLabels',[],'yLabel',[],'YTickLabels',[],'fontsize',16)
    else
        set(gca,'XTickLabels',[],'yLabel',[],'fontsize',16)
    end
end
saveas(gcf,[figures_dir 'allExp_' tag],'fig')
saveas(gcf,[figures_dir 'allExp_' tag],'pdf')
saveas(gcf,[figures_dir 'allExp_' tag],'png')

%% get subject lists
load([results_dir 'Tb'])
clear Subjects UIDs
for iexp=1:length(expNames)
    Subjects{iexp} = unique(Tb(ismember(Tb.ExpName,expNames{iexp}),:).Subject);
    UIDs{iexp} = unique(Tb(ismember(Tb.ExpName,expNames{iexp}),:).UID);
    disp(expNs{iexp})
    disp(sort(UIDs{iexp}))
end

