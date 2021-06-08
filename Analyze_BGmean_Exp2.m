%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir = [root_dir 'figures/'];
methods_dir = [root_dir 'Methods/'];

expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'};
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')
WordLikeness = [1.25,1.75,2.5,3.75,4.75];
%% all together

T=readtable([methods_dir 'Exp2/exp2_BGmean.xlsx']);
cols=7:3:40;
allN = nan(360,5);
for cond=1:5
    N = T{T.CondNum==cond,cols};N=N(:);
    allN(:,cond) = N;
    meanN(cond) = mean(N);
    steN(cond) = std(N)/sqrt(length(N));
end
BGmean=meanN;
save([methods_dir 'Exp2/exp2_BGmean'],'BGmean')


figure
%errorbar(meanN(end:-1:1),steN(end:-1:1),'linewidth',2)

errorbar(WordLikeness,meanN(end:-1:1),steN(end:-1:1),'linewidth',2)
%set(gca,'xtick',[5:-1:1],'xticklabel',{'Most','','','','Least'},'fontsize',16)
set(gca,'xtick',WordLikeness,'xticklabel',WordLikeness(end:-1:1),'fontsize',16)
set(gca,'ytick',[3000,4000],'fontsize',16)
ylim([2500 4500])

xlabel('Word-likeness (rating bin centers)')
ylabel('BG mean')
title('Phonotactics as a function of word-likeness')

saveas(gcf,[figures_dir 'Exp2_BGmean.fig'])
saveas(gcf,[figures_dir 'Exp2_BGmean.png'])
saveas(gcf,[figures_dir 'Exp2_BGmean.pdf'])

WordLikeness = WordLikeness';
meanN = meanN';
WordLikenessGroup = [1:5]';
T=table;
for cond=1:5
    for item = 1:size(allN,1)
        Temp = table;
        Temp.bgMean = allN(item,cond);
        Temp.WordLikeness = WordLikeness(cond);
        Temp.WordLikenessGroup = WordLikenessGroup(cond);
        T=[T;Temp];
    end
end
lm=fitlme(T,'bgMean ~ WordLikeness') 
lm.anova
%lme=fitlme(T,'bgMean ~ WordLikeness');
%cm=fitlme(T,'bgMean ~ 1');
%comp=compare(cm,lm)

[R,P]=corrcoef(T.WordLikeness,T.bgMean);

%% compare BGmean of two lowest conditions

disp('compare two least word-like groups')
[H,P,CI,STATS]  = ttest2(allN(:,1),allN(:,2))