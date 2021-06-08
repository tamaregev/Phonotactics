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
load([methods_dir 'Exp2/exp2_OrthoN'])
load([methods_dir 'Exp2/exp2_BGmean'])
[r, p]=corrcoef(OrthoN,BGmean);

%% prep table 

Tec=readtable([results_dir 'Exp2/Results.xls']);
Tecr = Tec(ismember(Tec.ROI,categorical([1:5])),:);

EffectName = Tecr.Effect;
Effect = zeros(size(EffectName));
Effect(strcmp(EffectName,'W')) = 2;
Effect(strcmp(EffectName,'W1')) = 1;
Effect(strcmp(EffectName,'W2')) = 0;
Effect(strcmp(EffectName,'W3')) = -1;
Effect(strcmp(EffectName,'W4')) = -2;
Tecr.EffectName = Tecr.Effect;
Tecr.Effect = Effect;

EffectName = Tecr.EffectName;
orthoN = zeros(size(EffectName));
orthoN(strcmp(EffectName,'W')) = OrthoN(5);
orthoN(strcmp(EffectName,'W1')) = OrthoN(4);
orthoN(strcmp(EffectName,'W2')) = OrthoN(3);
orthoN(strcmp(EffectName,'W3')) = OrthoN(2);
orthoN(strcmp(EffectName,'W4')) = OrthoN(1);
Tecr.OrthoN = orthoN;

EffectName = Tecr.EffectName;
bg = zeros(size(EffectName));
bg(strcmp(EffectName,'W')) = BGmean(5);
bg(strcmp(EffectName,'W1')) = BGmean(4);
bg(strcmp(EffectName,'W2')) = BGmean(3);
bg(strcmp(EffectName,'W3')) = BGmean(2);
bg(strcmp(EffectName,'W4')) = BGmean(1);
Tecr.bgMean = bg;


%% LME
formulas = {'EffectSize ~ bgMean + (1|ROI)'...
            'EffectSize ~ OrthoN + bgMean + (1|ROI)'...
            'EffectSize ~ OrthoN*bgMean  + (1|ROI)'...
            'EffectSize ~ OrthoN*bgMean + (OrthoN*bgMean|ROI)'...
};

lmes = cell(size(formulas));
for ilme=1:length(lmes)
    disp(ilme)
    lmes{ilme} = fitlme(Tecr,formulas{ilme});
    disp(lmes{ilme})
end

%% zscore measures

Tz = Tecr;
Tz.EffectSize = zscore(Tz.EffectSize);
Tz.OrthoN = zscore(Tz.OrthoN);
Tz.bgMean = zscore(Tz.bgMean);

lmesz = cell(size(formulas));
for ilme=1:length(lmesz)
    disp(ilme)
    lmesz{ilme} = fitlme(Tz,formulas{ilme},'Verbose',true);
    disp(lmesz{ilme})
end


%% test for main effects (Levy)
form = 'EffectSize ~ OrthoN*bgMean + (1 | ROI)';
lme = fitlme(Tz,form);

formNoBG = 'EffectSize ~ OrthoN*bgMean - bgMean + (1|ROI)';
formNoOrtho = 'EffectSize ~ OrthoN*bgMean - OrthoN + (1|ROI)';
lmeNoBG = fitlme(Tz,formNoBG);
lmeNoOrtho = fitlme(Tz,formNoOrtho);

disp('No BGmean')
compare(lmeNoBG,lme)

disp('No OrthoN')
compare(lmeNoOrtho,lme)

% for the full model with random terms:
formfull = 'EffectSize ~ OrthoN*bgMean + (OrthoN*bgMean | ROI)';
lmefull = fitlme(Tz,formfull);

formNoBGfull = 'EffectSize ~ OrthoN*bgMean - bgMean + (OrthoN*bgMean - bgMean |ROI)';
formNoOrthofull = 'EffectSize ~ OrthoN*bgMean - OrthoN + (OrthoN*bgMean - OrthoN|ROI)';

lmeNoBGfull = fitlme(Tz,formNoBGfull);
lmeNoOrthofull = fitlme(Tz,formNoOrthofull);

disp('No BGmean full')
compare(lmeNoBGfull,lmefull)

disp('No OrthoN full')
compare(lmeNoOrthofull,lmefull)

