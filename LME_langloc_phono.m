%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir   = [root_dir 'figures/'];
expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'};
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')
ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
conditions = { 'S','N'};
  
addpath('/Users/tamaregev/Dropbox/MATLAB/functions/myFunctions/')
filename = 'spm_ss_mROI_data.csv';

%% Read data
clear T
for iexp=1:length(expNames)
    expName = expNames{iexp};
    expN = expNs{iexp};
    T{iexp}=readtable([results_dir 'SN_Activations_Phonology' filesep expName '_SN' filesep filename]);
    switch expName
        case 'langlocSN'
            T{iexp} = T{iexp}(:,[3,2,4,5]);
        otherwise
            T{iexp} = T{iexp}(:,[1:3,5]);
    end
    ExpNum=iexp*ones(height(T{iexp}),1);
    ExpName = cell(size(ExpNum));
    ExpName(:,:) = {expName}; 
    T{iexp}=[T{iexp} table(ExpNum) table(ExpName)];
       UID=nan(height(T{iexp}),1);
       for ii=1:height(T{iexp})
           UID(ii)=str2num(T{iexp}.Subject{ii}(1:3));
       end
       T{iexp}=[table(UID), T{iexp}];
    ExpN = cell(size(ExpNum));
    ExpN(:,:) = {expN}; 
    T{iexp}=[T{iexp} table(ExpN)];
end
Tb = [T{1}; T{2}; T{3}; T{4}; T{5}];
save([results_dir 'Tb_langloc'],'Tb')

%%  each fROI
Tr = cell(length(expNames),5);%exp, roi

Tres=table;%a table that aggregates all results
for iexp=1:length(expNames)
    expName = expNames{iexp};
    disp(expName)
    expN = expNs{iexp};
    conds = conditions;
    
    Te=T{iexp};Te.UID=categorical(Te.UID);Te.ROI=categorical(Te.ROI);
%    effects = unique(Te.Effect);

  Tec = Te;
       EffectName = Tec.Effect;
       Effect = zeros(size(EffectName));
       Effect(strcmp(EffectName,'S')) = 0.5;
       Effect(strcmp(EffectName,'N')) = -0.5;
       Tec.Effect = Effect;
    formula = 'EffectSize ~ Effect + (1|UID)';
    for iroi = 1:5
        %disp([iroi ROIstring{iroi}])
        Tecr = Tec(ismember(Tec.ROI,categorical(iroi)),:);
        lme = fitlme(Tecr,formula);
        Tr{iexp,iroi} = lme2table(lme);

        P(iroi) = lme.anova.pValue(2);
        %anova(lme)
    end
    [P_FDR, H] = FDR(P,0.05);
    for iroi = 1:5
        disp([iroi ROIstring{iroi}])
        disp(['P=' num2str(P(iroi)) ', P_FDR=' num2str(P_FDR(iroi)) ' H=' num2str(H(iroi))])
        %save all results to a table
        Experiment = repmat(expNames(iexp),[size(Tr{iexp,iroi},1),1]);
        ROI = repmat(ROIstring(iroi),[size(Tr{iexp,iroi},1),1]);
        p_FDR = nan(size(Tr{iexp,iroi},1),1);p_FDR(2) = P_FDR(iroi);
        H_FDR = nan(size(Tr{iexp,iroi},1),1);H_FDR(2) = H(iroi);
        clear Tre
        Tre = [table(Experiment),table(ROI),Tr{iexp,iroi}(:,1:6),table(p_FDR),table(H_FDR),Tr{iexp,iroi}(:,7:end)];
        Tre.Properties.VariableNames{3} = 'FixedEffectName';
        Tres = [Tres; Tre];
    end
end

writetable(Tres,[results_dir filesep 'lmeLangloc_perROI.csv'])
