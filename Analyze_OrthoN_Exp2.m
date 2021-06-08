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
T=readtable([methods_dir 'Exp2/exp2_orthoN.xlsx']);
cols=6:2:28;% T(1,cols)
allN = nan(360,5);% n stimuli per condition x number of conditions
for cond=1:5
    N = T{T.CondNum==cond,cols};N=N(:);
    allN(:,cond) = N;
    meanN(cond) = mean(N);
    steN(cond) = std(N)/sqrt(length(N));
end

OrthoN=meanN;
save([methods_dir 'Exp2/exp2_OrthoN'],'OrthoN')

figure
errorbar(WordLikeness,meanN(end:-1:1),steN(end:-1:1),'linewidth',2)
set(gca,'xtick',WordLikeness,'xticklabel',WordLikeness(end:-1:1),'fontsize',16)
set(gca,'ytick',[0,0.4,0.8],'yticklabel',{'0','0.4','0.8'},'fontsize',16)

xlabel('Word-likeness')
ylabel('Orthographic neighborhood')
title('Exp 2')
saveas(gcf,[figures_dir 'Exp2_OrthoN.fig'])
saveas(gcf,[figures_dir 'Exp2_OrthoN.png'])

saveas(gcf,[figures_dir 'Exp2_OrthoN.pdf'])

disp('compare two least word-like groups')
[H,P,CI,STATS]  = ttest2(allN(:,1),allN(:,2))