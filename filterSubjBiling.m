%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir = [root_dir 'figures/'];
expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'};
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')
addpath('/Users/tamaregev/Dropbox/MATLAB/functions/myFunctions')

%% Read and prep data
iexp=3;%Biling
expName = expNames{iexp};
expN = expNs{iexp};
Te=readtable([results_dir expN filesep expName '_S.csv']);

UID=nan(height(Te),1);
for ii=1:height(Te)
   UID(ii)=str2num(Te.Subject{ii}(1:3));
end
Te=[table(UID), Te];
%%

iexp=3;
expName = expNames{iexp};
disp(expName)
expN = expNs{iexp};
conds = 'ENGW';
    
Te.UID=categorical(Te.UID);Te.ROI=categorical(Te.ROI);
Tec=Te(strcmp(Te.Effect,conds),:);
Tecr = Tec(ismember(Tec.ROI,categorical([1:5])),:);
% 
% formula = 'EffectSize ~ UID - 1 + (1|ROI)';
% lme = fitlme(Tecr,formula);
Ts = cell(numel(unique(Te.UID)),1);
for is = 1:numel(unique(Te.UID))
    UIDs = unique(Te.UID);
    UID = UIDs(is);disp(UID)
    Tecrs = Tecr(ismember(Tecr.UID,categorical(UID)),:);
    means(is)=mean(Tecrs.EffectSize);
    disp(['mean = ' num2str(means(is))])
    formula = 'EffectSize ~ 1 + (1|ROI)';
    lmes = fitlme(Tecrs,formula);
    Ts{is} = lme2table(lmes);
    
    disp(Ts{is})
end
%% select subjects
isNeg = nan(numel(unique(Te.UID)),1);
isNegSig = nan(numel(unique(Te.UID)),1);
for is = 1:numel(unique(Te.UID))
    Ttemp = Ts{is};
    isNeg(is) = Ttemp.Estimate < 0;
    isNegSig(is) = Ttemp.Estimate < 0 && Ttemp.pValue < 0.05 ;
end
UIDs = unique(Te.UID);
negatives = UIDs(find(isNeg));
negSignificants = UIDs(find(isNegSig));
disp(['S estimate is negative: ' negatives'])
disp(['S estimate is significant: ' negSignificants'])



%% write results into a csv

Ta = table;
for is = 1:numel(unique(Te.UID))
    Ta = [Ta;Ts{is}]; 
end
UID = unique(Te.UID);
Ta = [table(UID),Ta,table(isNeg),table(isNegSig)];
writetable(Ta,[results_dir expN filesep 'Stest.csv'])

