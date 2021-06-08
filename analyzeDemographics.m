%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir = [root_dir 'figures/'];
expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'};
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')
Excluded = {[],[],{'389','456','534'},[],[]};
%% read subject data
Tsa=readtable([results_dir 'Regev_etal_phonology_SUBJECTS.xlsx']);
Tsa=Tsa(1:end-6,:);
Tsa.Properties.VariableNames{4} = 'Biling';
% exclude participants
Tsa.Biling(ismember(Tsa.UID,Excluded{3}))=0;
%% N
N = nan(5,1);
for iexp = 1:5
    N(iexp) = sum(Tsa.(expNames{iexp})==1);
end
%% overlapping
overlap = nan(5,5);
for iexp = 1:5
    for jexp = 1:iexp
        overlap(iexp,jexp) = sum(ismember(Tsa.UID(Tsa.(expNames{iexp})==1),Tsa.UID(Tsa.(expNames{jexp})==1)));
    end
end
%% Females
F = nan(5,1);
for iexp = 1:5
     F(iexp)=sum(strcmp(Tsa.Sex(Tsa.(expNames{iexp})==1),'F'));
end

%% age
meanAge = nan(5,1);
stdAge = nan(5,1);
for iexp = 1:5
     meanAge(iexp)=mean(Tsa.Age_floor(Tsa.(expNames{iexp})==1));
     stdAge(iexp)=std(Tsa.Age_floor(Tsa.(expNames{iexp})==1));    
end
%% Left handedness
LHand = nan(5,1);
Amb = nan(5,1);
for iexp = 1:5
     LHand(iexp)=sum(strcmp(Tsa.Handedness(Tsa.(expNames{iexp})==1),'Left'));
     Amb(iexp)=sum(strcmp(Tsa.Handedness(Tsa.(expNames{iexp})==1),'Amb'));
end

%analysis per unique participants:
uniqueIDs = unique(Tsa.UID(Tsa.langlocSN | Tsa.SWJNaud | Tsa.Biling | Tsa.ParamNew | Tsa.Nhood));
sum(strcmp(Tsa.Handedness(Tsa.langlocSN | Tsa.SWJNaud | Tsa.Biling | Tsa.ParamNew | Tsa.Nhood),'Left'))
