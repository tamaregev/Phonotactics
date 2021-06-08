%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir   = [root_dir 'figures/'];
expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'};
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')
ROIstring = {'IFGorb','IFG','MFG','AntTemp','PostTemp','AngG','IFGorb-R','IFG-R','MFG-R','AntTemp-R','PostTemp-R','AngG-R'};
conditions = { 'N langloc','N aud','ENGNW',{'W','W1','W2','W3','W4'},{'high','low'}};
addpath('/Users/tamaregev/Dropbox/MATLAB/functions/myFunctions/')
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

%% critical experiments
%% each fROI
Tr = cell(length(expNames),5);%exp, roi

Tres=table;% a table that aggregates all results
for iexp=1:length(expNames)
    expName = expNames{iexp};
    disp(expName)
    expN = expNs{iexp};
    conds = conditions{iexp};
    
    Te=T{iexp};Te.UID=categorical(Te.UID);Te.ROI=categorical(Te.ROI);
%    effects = unique(Te.Effect);
    if ~iscell(conds) % only 1 level for the effect
        Tec=Te(strcmp(Te.Effect,conds),:);
        formula = 'EffectSize ~ 1 + (1|UID)';

        for iroi = 1:5
            %disp([iroi ROIstring{iroi}]);
            Tecr = Tec(ismember(Tec.ROI,categorical(iroi)),:);
            lme = fitlme(Tecr,formula);
            Tr{iexp,iroi} = lme2table(lme);
            
            P(iroi) = lme.Coefficients.pValue;
            Es(iroi) = lme.Coefficients.Estimate;
            
        end
        
        [P_FDR, H] = FDR(P,0.05);
        
        for iroi = 1:5
            disp([iroi ROIstring{iroi}])
            disp(['Estimate=' num2str(Es(iroi)) ', P=' num2str(P(iroi)) ', P_FDR=' num2str(P_FDR(iroi)) ' H=' num2str(H(iroi))])
            %save all results to a table
            Experiment = repmat(expNames(iexp),[size(Tr{iexp,iroi},1),1]);
            ROI = repmat(ROIstring(iroi),[size(Tr{iexp,iroi},1),1]);
            p_FDR = repmat(P_FDR(iroi),[size(Tr{iexp,iroi},1),1]);
            H_FDR = repmat(H(iroi),[size(Tr{iexp,iroi},1),1]);
            clear Tre
            Tre = [table(Experiment),table(ROI),Tr{iexp,iroi}(:,1:6),table(p_FDR),table(H_FDR),Tr{iexp,iroi}(:,7:end)];
            Tre.Properties.VariableNames{3} = 'FixedEffectName';
            Tres = [Tres; Tre];
        end

        
    else % effect has several levels
        logical = zeros(size(Te.Effect));
        for ic=1:length(conds)
            logical = logical | strcmp(Te.Effect,conds{ic});
            %disp(sum(logical))
        end
        Tec=Te(logical,:);
        if iexp==4
           EffectName = Tec.Effect;
           Effect = zeros(size(EffectName));
           Effect(strcmp(EffectName,'W')) = 2;
           Effect(strcmp(EffectName,'W1')) = 1;
           Effect(strcmp(EffectName,'W2')) = 0;
           Effect(strcmp(EffectName,'W3')) = -1;
           Effect(strcmp(EffectName,'W4')) = -2;
           Tec.Effect = Effect;
        elseif iexp==5
           EffectName = Tec.Effect;
           Effect = zeros(size(EffectName));
           Effect(strcmp(EffectName,'low')) = -0.5;
           Effect(strcmp(EffectName,'high')) = 0.5;
           Tec.Effect = Effect;
        end
        formula = 'EffectSize ~ Effect + (Effect|UID)';
        for iroi = 1:5
            %disp([iroi ROIstring{iroi}])
            Tecr = Tec(ismember(Tec.ROI,categorical(iroi)),:);
            lme = fitlme(Tecr,formula);
            Tr{iexp,iroi} = lme2table(lme);

            P_INT(iroi) = lme.anova.pValue(1);
            
            P_SLOPE(iroi) = lme.anova.pValue(2);
            %anova(lme)
        end
        [P_INT_FDR, H_INT] = FDR(P_INT,0.05);
        [P_SLOPE_FDR, H_SLOPE] = FDR(P_SLOPE,0.05);

        for iroi = 1:5
            disp([iroi ROIstring{iroi}])
            disp(['P_INT=' num2str(P_INT(iroi)) ', P_INT_FDR=' num2str(P_INT_FDR(iroi)) ' H=' num2str(H_INT(iroi))])
            disp(['P_SLOPE=' num2str(P_INT(iroi)) ', P_SLOPE_FDR=' num2str(P_SLOPE_FDR(iroi)) ' H=' num2str(H_SLOPE(iroi))])

            %save all results to a table
            Experiment = repmat(expNames(iexp),[size(Tr{iexp,iroi},1),1]);
            ROI = repmat(ROIstring(iroi),[size(Tr{iexp,iroi},1),1]);
            p_FDR = nan(size(Tr{iexp,iroi},1),1);p_FDR(1) = P_INT_FDR(iroi);p_FDR(2) = P_SLOPE_FDR(iroi);
            H_FDR = nan(size(Tr{iexp,iroi},1),1);H_FDR(1) = H_INT(iroi);H_FDR(2) = H_SLOPE(iroi);
            clear Tre
            Tre = [table(Experiment),table(ROI),Tr{iexp,iroi}(:,1:6),table(p_FDR),table(H_FDR),Tr{iexp,iroi}(:,7:end)];
            Tre.Properties.VariableNames{3} = 'FixedEffectName';
            Tres = [Tres; Tre];
        end
    end
end

writetable(Tres,[results_dir filesep 'lmeCrit_perROI.csv'])

%%  all fROIs together
Tr = cell(length(expNames));%exp
Tres=table;%a table that aggregates all results

for iexp=1:length(expNames)
    expName = expNames{iexp};
    disp(expName)
    expN = expNs{iexp};
    conds = conditions{iexp};
    
    Te=T{iexp};Te.UID=categorical(Te.UID);Te.ROI=categorical(Te.ROI);
%    effects = unique(Te.Effect);
    if ~iscell(conds)%effect has only 1 level = intercept
        Tec=Te(strcmp(Te.Effect,conds),:);
        Tecr = Tec(ismember(Tec.ROI,categorical([1:5])),:);
        formula = 'EffectSize ~ 1 + (1|UID) + (1|ROI)';

        lme = fitlme(Tecr,formula);
        Tr{iexp} = lme2table(lme);

        anova(lme);
        %    [H,P,CI,STATS] = ttest(Tecr.EffectSize);
        %    disp(['t=(' num2str(STATS.df) ')' num2str(STATS.tstat) ', p=' num2str(P)])
        P = lme.Coefficients.pValue;
        es = lme.Coefficients.Estimate;
        disp(['Estimate=' num2str(es) , 'P=' num2str(P) ' H=' num2str(P<0.05)])
            %save all results to a table
            Experiment = repmat(expNames(iexp),[size(Tr{iexp},1),1]);
            H = P<0.05;
            clear Tre
            Tre = [table(Experiment),Tr{iexp}(:,1:6),table(H),Tr{iexp}(:,7:end)];
            Tre.Properties.VariableNames{2} = 'FixedEffectName';
            Tres = [Tres; Tre];
            
    else %effect has several levels
        logical = zeros(size(Te.Effect));
        for ic=1:length(conds)
            logical = logical | strcmp(Te.Effect,conds{ic});
            %disp(sum(logical))
        end
        Tec=Te(logical,:);
        if iexp==4
           EffectName = Tec.Effect;
           Effect = zeros(size(EffectName));
%            Effect(strcmp(EffectName,'W')) = -2;
%            Effect(strcmp(EffectName,'W1')) = -1;
%            Effect(strcmp(EffectName,'W2')) = 0;
%            Effect(strcmp(EffectName,'W3')) = 1;
%            Effect(strcmp(EffectName,'W4')) = 2;
           Effect(strcmp(EffectName,'W')) = 2;
           Effect(strcmp(EffectName,'W1')) = 1;
           Effect(strcmp(EffectName,'W2')) = 0;
           Effect(strcmp(EffectName,'W3')) = -1;
           Effect(strcmp(EffectName,'W4')) = -2;
           Tec.Effect = Effect;
        elseif iexp==5
           EffectName = Tec.Effect;
           Effect = zeros(size(EffectName));
           Effect(strcmp(EffectName,'low')) = -0.5;
           Effect(strcmp(EffectName,'high')) = 0.5;
           Tec.Effect = Effect;
        end
        Tecr = Tec(ismember(Tec.ROI,categorical([1:5])),:);
        formula = 'EffectSize ~ Effect + (Effect|UID) + (Effect|ROI)';
        lme = fitlme(Tecr,formula);
        Tr{iexp} = lme2table(lme);

        anova(lme);
        P_INT = lme.anova.pValue(1);
        P_SLOPE = lme.anova.pValue(2);
        
        %es = lme.Coefficients.Estimate;
        disp(['P_INT=' num2str(P_INT) ' H_INT=' num2str(P_INT<0.05)])
        disp(['P_SLOPE=' num2str(P_SLOPE) ' H_SLOPE=' num2str(P_SLOPE<0.05)])
       
        %save all results to a table
            Experiment = repmat(expNames(iexp),[size(Tr{iexp},1),1]);
            H = nan(size(Tr{iexp},1),1);H(1) = P_INT<0.05;H(2) = P_SLOPE<0.05;
            clear Tre
            Tre = [table(Experiment),Tr{iexp}(:,1:6),table(H),Tr{iexp}(:,7:end)];
            Tre.Properties.VariableNames{2} = 'FixedEffectName';
            Tres = [Tres; Tre];
    end
    

end
writetable(Tres,[results_dir filesep 'lmeCrit.csv'])

%% calculating for the auditory experiments together
   
%
disp('auditory experiments (1b and 1c) together')
Te2=T{2};conds2 = conditions{2};    
Tec2=Te2(strcmp(Te2.Effect,conds2),:);
Te3=T{3};conds3 = conditions{3};
Tec3=Te3(strcmp(Te3.Effect,conds3),:);
    
Tec=[Tec2;Tec3];
Tec.UID=categorical(Tec.UID);Tec.ROI=categorical(Tec.ROI);

    
formula = 'EffectSize ~ 1 + (1|UID) + (1|ROI) + (1|Effect)';%Effect here encodes the experiment
%all ROIs together
disp('all ROIs')
lme = fitlme(Tec,formula);
anova(lme);
P = lme.Coefficients.pValue;
es = lme.Coefficients.Estimate;
disp(['Estimate=' num2str(es) , 'P=' num2str(P) ' H=' num2str(P<0.05)])


%each ROI separately
disp('individual ROIs')

formula = 'EffectSize ~ 1 + (1|UID) + (1|Effect)';

for iroi = 1:5
%    disp([iroi ROIstring{iroi}])
    Tecr = Tec(ismember(Tec.ROI,categorical(iroi)),:);
    lme = fitlme(Tecr,formula);
    anova(lme);
    P(iroi) = lme.Coefficients.pValue;
    Es(iroi) = lme.Coefficients.Estimate;
%    [H,P,CI,STATS] = ttest(Tecr.EffectSize);
%    disp(['t=(' num2str(STATS.df) ')' num2str(STATS.tstat) ', p=' num2str(P)])
end  

[P_FDR, H] = FDR(P,0.05);
for iroi = 1:5
    disp([iroi ROIstring{iroi}])
    disp(['Estimate = ' num2str(Es(iroi)) ', P=' num2str(P(iroi)) ', P_FDR=' num2str(P_FDR(iroi)) ' H=' num2str(H(iroi))])
end


