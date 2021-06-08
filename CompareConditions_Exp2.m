%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir = [root_dir 'figures/'];
methods_dir = [root_dir 'Methods/'];

expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'};
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')
WordLikeness = [1.25,1.75,2.5,3.75,4.75];
%%
Tec=readtable([results_dir 'Exp2/Results.xls']);
Tecr = Tec(ismember(Tec.ROI,categorical([1:5])),:);

if 0 % reduce to include only W3 and W4
%     EffectName = Tecr.Effect;
%     Effect = zeros(size(EffectName));
%     Effect(strcmp(EffectName,'W')) = 5;
%     Effect(strcmp(EffectName,'W1')) = 4;
%     Effect(strcmp(EffectName,'W2')) = 3;
%     Effect(strcmp(EffectName,'W3')) = 1;
%     Effect(strcmp(EffectName,'W4')) = -1;
%     Tecr.Effect = Effect;
%     T = Tecr(ismember(Tecr.Effect,[1 -1]),:);
W34 = strcmp(Tecr.Effect,{'W3'}) | strcmp(Tecr.Effect,{'W4'});
T = Tecr(W34,:);
else
    T=Tecr;
    T.Effect = categorical(T.Effect);
    C = categories(T.Effect);
    C2 = C(end:-1:1);
    T.Effect = reordercats(T.Effect,C2);
end
formula = 'EffectSize ~ Effect + (1|UID) + (1|ROI)';
lme = fitlme(T,formula)

%% Comparing to other models:
formula = 'EffectSize ~ 1 + (1|UID) + (1|ROI)';
lme_const = fitlme(T,formula,'Verbose',true);
compare(lme_const,lme)

%the more complex models are not significantly different:
formula = 'EffectSize ~ Effect + (Effect|UID) + (1|ROI)';
lme_UIDslope = fitlme(T,formula,'Verbose',true);

compare(lme,lme_UIDslope)

formula = 'EffectSize ~ Effect + (Effect|UID) + (Effect|ROI)';
lme_ROIslope = fitlme(T,formula,'Verbose',true)

compare(lme_UIDslope,lme_ROIslope)

%% coding as a slope var:
Tec=readtable([results_dir 'Exp2/Results.xls']);
Tecr = Tec(ismember(Tec.ROI,categorical([1:5])),:);
W34 = strcmp(Tecr.Effect,{'W3'}) | strcmp(Tecr.Effect,{'W4'});
T = Tecr(W34,:);

     EffectName = T.Effect;
     EffectNum = zeros(size(EffectName));

     EffectNum(strcmp(EffectName,'W3')) = 0.5;
     EffectNum(strcmp(EffectName,'W4')) = -0.5;
     T.EffectNum = EffectNum;
    

formula = 'EffectSize ~ EffectNum + (1|UID) + (1|ROI)';
lme = fitlme(T,formula)
