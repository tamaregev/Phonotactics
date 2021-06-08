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
T=readtable([methods_dir 'Exp3/exp3_orthoN.xls']);

cols=[6,17];
clear meanN steN
for cond=1:2
    allN = T{:,cols(cond)};
    meanN(cond) = mean(allN);
    steN(cond) = std(allN)/sqrt(length(allN));
end
figure
errorbar(meanN,steN,'linewidth',2,'Color','m')
xlim([0.75 2.25])
set(gca,'xtick',[1:2],'XtickLabel',{'high','low'},'fontsize',16)
set(gca,'ytick',[1:5],'fontsize',16)

xlabel('Nhood')
ylabel('Orthographic neighborhood')
title('Exp 3')
saveas(gcf,['/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/figures/Exp3_OrthoN.fig'])
saveas(gcf,['/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/figures/Exp3_OrthoN.png'])
saveas(gcf,['/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/figures/Exp3_OrthoN.pdf'])
[H,P,CI,STATS]  = ttest2(T{:,cols(1)},T{:,cols(2)})
%%
[R,P] = corrcoef([T.neighbors_x ; T.neighbors_y ],[T.ortho_nx; T.ortho_ny]);
