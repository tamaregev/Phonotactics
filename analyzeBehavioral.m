%% Definitions
root_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/';
results_dir = [root_dir 'Results/'];
figures_dir = [root_dir 'figures/'];
expNames = {'langlocSN','SWJNaud','Biling','ParamNew','Nhood'};
expNs = {'Exp1a','Exp1b','Exp1c','Exp2','Exp3'}; 
beh_filename = 'Regev_etal_phonology_BEHAVIORAL.xlsx';
addpath(genpath('/Users/tamaregev/Dropbox/MATLAB/functions/downloaded'));

%% langloc_SN
iexp = 1;
%load stim order:
load('/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/Methods/Exp1a/mainParvizi_langloc2013-2017_SN2cond_stimSets/langloc_fmri_run1_stim_v1.mat')
order{1}=stim(:,end);
load('/Users/tamaregev/Dropbox/postdoc/Fedorenko/Phonotactics/Methods/Exp1a/mainParvizi_langloc2013-2017_SN2cond_stimSets/langloc_fmri_run2_stim_v1.mat')
order{2}=stim(:,end);

folder = [results_dir expNs{iexp} filesep 'behavioral' filesep 'Parvizi_behavioral'];
D = dir(folder);
dirFlags = [D.isdir];
D = D(dirFlags);
D = D(~ismember({D.name}, {'.', '..','.DS_Store'})); % dir returns '.' and '..' (usually in first slot)
Ta = table;
for id=1:numel(D)
    Di = dir([folder filesep D(id).name]);
    Di = Di(~ismember({Di.name}, {'.', '..','.DS_Store'})); % dir returns '.' and '..' (usually in first slot)
    S_did=0;
    N_did=0;
    S_rt=nan(96,1);%j=0;
    N_rt=nan(96,1);%m=0;
    isset = 'nn';
    reversed = 'nn';
    runs = 'nn';
    
    for iid=1:numel(Di)% the 2 runs
        Temp = table;
        load([Di(iid).folder filesep Di(iid).name])%loads subj_data
        runn=str2double(Di(iid).name(strfind(Di(iid).name,'run')+3));
        if strfind(Di(iid).name,'set'), isset(iid) = Di(iid).name(strfind(Di(iid).name,'set')+3); else isset(iid) = 'v'; end
        if isfield(subj_data,'reversed')
            if subj_data.reversed, reversed(iid)='1'; else, reversed(iid)='0'; end
        end
        runs(iid) = num2str(subj_data.run);

        ord = order{runn};
        for itr=1:numel(subj_data.did_respond)
            switch ord{itr}
                case 'S'
                   if subj_data.did_respond(itr)==1
                       S_did = S_did+1;
                       S_rt(itr+(48*(iid-1))) = subj_data.rt(itr);
                   end
                case 'N'
                    if subj_data.did_respond(itr)==1
                        N_did = N_did+1;
                        N_rt(itr+(48*(iid-1))) = subj_data.rt(itr);                   
                    end
            end
        end
    end
    Temp.subj = {D(id).name};
    Temp.Nruns = numel(Di);
    Temp.runs = runs;
    Temp.sets = isset;%v if v1/v2 materials. otherwise gives the version of set# out of the 5
    Temp.reversed = reversed;

    Temp.S_did = S_did;
    Temp.N_did = N_did;
    Temp.S_rt = nanmean(S_rt);
    Temp.N_rt = nanmean(N_rt);
    Temp.S_rt_sd = nanstd(S_rt);
    Temp.N_rt_sd = nanstd(N_rt);
    Ta = [Ta; Temp];
end
writetable(Ta, [folder filesep 'Parvizi_behavioral.xls'])
S_did_mean = nanmean(Ta.S_did/48*100)
S_did_std = nanstd(Ta.S_did/48*100)

N_did_mean = nanmean(Ta.N_did/48*100)
N_did_std = nanstd(Ta.N_did/48*100)

[H,P,CI,STATS]  = ttest(Ta.S_did/48*100,Ta.N_did/48*100)

S_rt_mean = nanmean(Ta.S_rt*1000)
S_rt_std = nanstd(Ta.S_rt*1000)

N_rt_mean = nanmean(Ta.N_rt*1000)
N_rt_std = nanstd(Ta.N_rt*1000)

[H,P,CI,STATS]  = ttest(Ta.S_rt*1000,Ta.N_rt*1000)

% calcs for by block analysis:
subj_v=Ta.subj(strcmp(cellstr(Ta.sets),'vv'));
uniqueIDs = subj_v;

%% SWJNaud
iexp = 2;
conditions = {'S','W','J','N'};

%folder = [results_dir expNs{iexp} filesep 'behavioral'];
%T = readtable([folder filesep 'SWJN_audio.xls']);
%newer dataset:
folder = [results_dir expNs{iexp} filesep 'behavioral' filesep 'indiv-subj-folders'];
T = readtable([folder filesep 'SWJN_audio_new.xls']);


subjects = unique(T.Subj_ID);
RT = nan(numel(conditions),numel(subjects));
acc = nan(numel(conditions),numel(subjects));
for con = 1:numel(conditions)
    for is=1:numel(subjects)    
        Ntrials = height(T(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con ,:));
        % P = probe R = response, app = appeared noapp = didn't appear
        Papp_Rapp = height(T(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==1 & T.Response==1,:));
        Pnoapp_Rnoapp = height(T(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==0 & T.Response==2,:));
        Papp_Rnoapp = height(T(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==1 & T.Response==2 ,:));
        Pnoapp_Rapp = height(T(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==0 & T.Response==1,:));
        Papp_Rnull = height(T(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==1 & T.Response==0,:));
        Pnoapp_Rnull = height(T(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==0 & T.Response==0,:));
        
        RT_Papp_Rapp = T.RT(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==1 & T.Response==1,:);
        RT_Pnoapp_Rnoapp = T.RT(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==0 & T.Response==2,:);
        
%        acc(con,is) = correct / Ntrials;
        acc(con,is) = (Papp_Rapp + Pnoapp_Rnoapp) / Ntrials;
%        RT(con,is) = mean(T.RT(strcmp(T.Subj_ID,subjects{is}) & T.Condition==con & T.Mem_trial==1,:));
        RT(con,is) = (mean(RT_Papp_Rapp)+mean(RT_Pnoapp_Rnoapp))./2;
    end
end

%% ParamNew
iexp = 4;
WordLikeness = [1.25,1.75,2.5,3.75,4.75];

sheet = 'Exp2';
conditions = {'Word-list','Nonword-list_a','Nonword-list_b','Nonword-list_c','Nonword-list_d'};
T = readtable([results_dir beh_filename], 'Sheet',expNs{iexp});
% I ran this only once to eliminate irrelevant conditions from Mollica etal
%select just the conditions:
% logical = zeros(size(T.CondName));
% for ic=1:length(conditions)
%     logical = logical | strcmp(T.CondName,conditions{ic});
%     %disp(sum(logical))
% end
% T=T(logical,:);
% writetable(T,[results_dir beh_filename], 'Sheet',expNs{iexp})

subjects = unique(T.SubjectID);
Tr = table;

for is=1:length(subjects)
    for ic=1:length(conditions)
        for ip=1:2
            Temp = table;
            Temp.Subject = subjects(is);
            Temp.Condition = conditions(ic);
            Temp.ProbeAppeared = ip-1;
            Temp.N = sum(strcmp(T.SubjectID,subjects(is)) & strcmp(T.CondName,conditions(ic)) & ismember(T.ProbeAppeared,ip-1));
            Temp.Ncorrect = sum(strcmp(T.SubjectID,subjects(is)) & strcmp(T.CondName,conditions(ic)) & ismember(T.ProbeAppeared,ip-1) & T.AnsweredCorrectly==1);
            Temp.Nincorrect = sum(strcmp(T.SubjectID,subjects(is)) & strcmp(T.CondName,conditions(ic)) & ismember(T.ProbeAppeared,ip-1) & T.AnsweredCorrectly==0);
            Temp.RT_correct = nanmean(T.RT(strcmp(T.SubjectID,subjects(is)) & strcmp(T.CondName,conditions(ic)) & ismember(T.ProbeAppeared,ip-1) & T.AnsweredCorrectly==1));
            Temp.RT_incorrect = nanmean(T.RT(strcmp(T.SubjectID,subjects(is)) & strcmp(T.CondName,conditions(ic)) & ismember(T.ProbeAppeared,ip-1) & T.AnsweredCorrectly==0));
            Tr = [Tr; Temp];
        end
    end
end
Tr.Condition = categorical(Tr.Condition,conditions,'Ordinal',true);
conditions = categorical(conditions,conditions,'Ordinal',true);

[s,I]=sort(Tr.Condition);
Tr = Tr(I,:);
%
Hits = nan(numel(subjects),numel(conditions));
RT_Hits = nan(numel(subjects),numel(conditions));
FAs = nan(numel(subjects),numel(conditions));
RT_FAs = nan(numel(subjects),numel(conditions));
Miss = nan(numel(subjects),numel(conditions));
RT_Miss = nan(numel(subjects),numel(conditions));
CRs = nan(numel(subjects),numel(conditions));
RT_CRs = nan(numel(subjects),numel(conditions));
Corrects = nan(numel(subjects),numel(conditions));
RT_corrects = nan(numel(subjects),numel(conditions));

dps = nan(numel(subjects),numel(conditions));
crits = nan(numel(subjects),numel(conditions));
for is = 1:numel(subjects)
    for ic=1:numel(conditions)
        %errors:
        FAs(is,ic) = Tr.Nincorrect(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==0)/Tr.N(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==0);
        Miss(is,ic) = Tr.Nincorrect(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==1)/Tr.N(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==1);
        RT_FAs(is,ic) = Tr.RT_incorrect(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==0)/Tr.N(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==0);
        RT_Miss(is,ic) = Tr.RT_incorrect(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==1)/Tr.N(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==1);
       
        %correct:
        NcorrApp = Tr.Ncorrect(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==1);
        Napp = Tr.N(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==1);
        Hits(is,ic) = NcorrApp/Napp;
        NcorrNoApp = Tr.Ncorrect(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==0);
        NnoApp = Tr.N(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==0);
        CRs(is,ic) = NcorrNoApp/NnoApp;
        Corrects(is,ic) = (NcorrApp+NcorrNoApp) / (Napp+NnoApp); 
        
        RT_Hits(is,ic) = Tr.RT_correct(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==1);
        RT_CRs(is,ic) = Tr.RT_correct(strcmp(Tr.Subject,subjects(is)) & ismember(Tr.Condition,conditions(ic)) & Tr.ProbeAppeared==0);
        RT_corrects(is,ic) = (RT_Hits(is,ic)*NcorrApp + RT_CRs(is,ic)*NcorrNoApp)/(NcorrApp+NcorrNoApp);
        
        [dp,c] = dprime_simple(Hits(is,ic),FAs(is,ic));
        dps(is,ic) = dp;
        crits(is,ic) = c;        
    end
end
mean(Corrects)
std(Corrects)
conditions
mean(RT_corrects)
std(RT_corrects)


maxhits = 0.93;
minFAs = 0.04;
maxdp = dprime_simple(maxhits,minFAs);
dps(dps==Inf)=maxdp;

% plot measure
name = 'Corrects';
measure = eval(name);
meanmeasure = nanmean(measure);
stdmeasure = nanstd(measure);

figure
%plot(WordLikeness(end:-1:1),meanmeasure(end:-1:1),'linewidth',2)
%hold on
errorbar(WordLikeness(end:-1:1),meanmeasure(end:-1:1),stdmeasure(end:-1:1)./sqrt(numel(subjects)),'linewidth',2)
%set(gca,'xTick',1:5)
%set(gca,'xTickLabel',strrep(cellstr(conditions),'_',' '),'fontsize',16)
title(strrep(name,'_',' '))
set(gca,'xtick',WordLikeness,'xticklabel',WordLikeness(end:-1:1),'fontsize',16)
set(gca,'ytick',[0,0.5,1],'yticklabel',{'0%','50%','100%'},'fontsize',16)
ylim([0,1])

xlabel('Word-likeness (rating bin centers)')
ylabel('Accuracy (%)')
title('Behavioral responses in Exp. 2')

saveas(gcf,[figures_dir 'Exp2_behavior.fig'])
saveas(gcf,[figures_dir 'Exp2_behavior.png'])
saveas(gcf,[figures_dir 'Exp2_behavior.pdf'])

WordLikeness = WordLikeness';
%meanN = meanN';
WordLikenessGroup = [1:5]';
T=table;
for cond=1:5
    for subj = 1:size(measure,1)
        Temp = table;
        Temp.measure = measure(subj,cond);
        Temp.WordLikeness = WordLikeness(cond);
        Temp.subj = subj;
        T=[T;Temp];
    end
end
T.subj=categorical(T.subj);
lme = fitlme(T,'measure ~ WordLikeness + (WordLikeness|subj)'); 

%% Nhood
% read subject data
%Tsa=readtable([results_dir 'Regev_etal_phonology_SUBJECTS.xlsx']);

iexp=5;
load([results_dir 'Tb'],'Tb')
subjects = unique(Tb.Subject(Tb.ExpNum==iexp));
folder = [results_dir expNs{iexp} filesep 'behavioral'];
%rename files starting with FED2016* to FED_2016*
% badfiles=dir([folder filesep 'FED2016*']);
% for ib=1:numel(badfiles)
%     movefile([folder filesep badfiles(ib).name],[folder filesep 'FED_' badfiles(ib).name(4:end)])
% end

accuracy = nan(numel(subjects),2);
RT = nan(numel(subjects),2);
nruns=nan(numel(subjects),1);
for is=1:length(subjects)
    session = subjects{is};
    behfile = [session(5:end-7)];
    files = dir([folder filesep behfile '*']);
    maxnruns = 3;
    acc = nan(maxnruns,2);
    rt = nan(maxnruns,2);
    for irun=1:maxnruns
        for xn=1:7
            filename = [behfile '_run' num2str(irun) '-x' num2str(xn) '_data.csv'];
            if exist([folder filesep filename],'file')
                Tsr=readtable([folder filesep filename]);
                if height(Tsr)
                    %disp([session ' run=' num2str(irun) ' x=' num2str(xn)])
                    acc(irun,1)=sum(Tsr.accuracy(Tsr.N==0)==1)/numel(Tsr.accuracy(Tsr.N==0));
                    acc(irun,2)=sum(Tsr.accuracy(Tsr.N==1)==1)/numel(Tsr.accuracy(Tsr.N==1));
                    
                    rt(irun,1)=nanmean(Tsr.RT(Tsr.N==0));
                    rt(irun,2)=nanmean(Tsr.RT(Tsr.N==1));
                    
                end
            end
        end        
    end
    nruns(is) = sum(~isnan(acc(:,1)));
    accuracy(is,:) = nanmean(acc);
    RT(is,:) = nanmean(rt);
end

[Hrt,Prt,CIrt,STATSrt]  = ttest(RT(:,1),RT(:,2))

[Hac,Pac,CIac,STATSac]  = ttest(accuracy(:,1),accuracy(:,2))

figure
errorbar(mean(accuracy),std(accuracy)/sqrt(size(accuracy,1)),'linewidth',2,'Color','m')
ylim([0 1])
%set(gca,'ytick',[3000,4000],'fontsize',16)

set(gca,'ytick',[0 0.5 1],'YtickLabel',{'0','50','100'},'fontsize',16)
set(gca,'xtick',[1:2],'XtickLabel',{'high','low'},'fontsize',16)
xlim([0.75 2.25])
xlabel('Neighborhood size')
ylabel('Accuracy (%)')
title('Exp 3 - Behavioral')

saveas(gcf,[figures_dir 'Exp3_behavior.fig'])
saveas(gcf,[figures_dir 'Exp3_behavior.png'])
saveas(gcf,[figures_dir 'Exp3_behavior.pdf'])