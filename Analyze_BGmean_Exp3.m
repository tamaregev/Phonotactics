%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir = [root_dir 'figures/'];
methods_dir = [root_dir 'Methods/'];

expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'};
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions/mine')

%%

T=readtable([methods_dir 'Exp3/exp3_BGmean.xls']);

cols=[8,20];
clear meanN steN
for cond=1:2
    allN = T{:,cols(cond)};
    meanN(cond) = mean(allN);
    steN(cond) = std(allN)/sqrt(length(allN));
end
figure
errorbar(meanN,steN,'linewidth',2,'Color','m')
ylim([2500 4500])
set(gca,'ytick',[3000,4000],'fontsize',16)

set(gca,'xtick',[1:2],'XtickLabel',{'high','low'},'fontsize',16)
xlim([0.75 2.25])
xlabel('Nhood')
ylabel('BGmean')
title('Exp 3 - BGmean')


saveas(gcf,[figures_dir 'Exp3_BGmean.fig'])
saveas(gcf,[figures_dir 'Exp3_BGmean.png'])
saveas(gcf,[figures_dir 'Exp3_BGmean.pdf'])

[H,P,CI,STATS]  = ttest2(T{:,cols(1)},T{:,cols(2)})
