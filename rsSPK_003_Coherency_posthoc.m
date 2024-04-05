%% Coherency post-processing
%-
%- Require fieldtrip and to have run "..._002_LFPpreproc.m"
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2022.12

%% average across PRE/POST and all drug condition to extract FQ of interest

clear
sessions = {'L082522a' 'L090122a' 'L081822b' 'L090822a' 'L092222a' 'H041422a' 'H051222a' 'H042122a' 'H052422a' 'H060322a'}
cond = {'DCZ' 'DCZ'  'Saline' 'Saline' 'CNO' 'DCZ' 'DCZ' 'Saline' 'Saline' 'CNO'};

path2go = 'I:\Havoc_REC\';

load([path2go 'rsSPK_elecreject.mat']) %- load the rejection file

ref = 'bipol'
measure2take = {'coh' 'coh_perm'} %- can be 'powcorr', 'coh' or 'wpli_debiased'

subdata=cell(1,2);
subdata_perm=cell(1,2);
for s = 1 : length(sessions)
    list = dir([path2go sessions{s} '\*corr_permfilt_' ref '.mat']);
    if ~isempty(list)
        load([path2go sessions{s} '\' list(1).name])

        %- extract average connectivity across contact for each FQ
        dumm_mean_both=[];
        dumm_perm_mean_both= [];
        for fq = 1 : length(conn.freq)

            eval(['coh = mean(conn(1).' measure2take{1} '(:,:,fq),3);']);
            eval(['coh_POST = mean(conn_post.' measure2take{1} '(:,:,fq),3);']);
            eval(['coh_perm = mean(conn(1).' measure2take{2} '(:,:,fq),3);']);
            eval(['coh_perm_POST = mean(conn_post.' measure2take{2} '(:,:,fq),3);']);

            %- remove the diagonal (by adding a matrix of 0s with NaN on the diagonal!)
            v = NaN(length(coh),1);
            D = diag(v);
            coh = coh + D;
            coh_POST = coh_POST + D;
            coh_perm = coh_perm + D;
            coh_perm_POST = coh_perm_POST + D;

            if strcmp(sessions{s}(1),'L')
                grp = [ones(size(coh,1)/3,1) ; 2*ones(size(coh,1)/3,1) ; 3*ones(size(coh,1)/3,1)]; %- 2 elec with same contact number
                mk = 1;
            else
                grp = [ones(size(coh,1)/2,1) ; 2*ones(size(coh,1)/2,1)]; %- 2 elec with same contact number
                mk = 2 ;
            end

            %- remove elc that were identified as too noisy
            sess_rej = find(ismember(session_reject,sessions{s})); %- match current session with rejection matrix
            keep_elec = true(length(coh),1);
            keep_elec(elec_reject{sess_rej})=false;

            % extract elec groups
            nbgrp = unique(grp);
            nbgrp(isnan(nbgrp))=[];

            for g1 = 1 : length(nbgrp)
                for g2 = 1 : length(nbgrp)
                    if g1>=g2 %- only do that below diagonal, to remove the symetrical matrix
                        du = coh(grp==g1 & keep_elec,grp==g2 & keep_elec);
                        du_post = coh_POST(grp==g1 & keep_elec,grp==g2 & keep_elec);
                        du_perm = coh_perm(grp==g1 & keep_elec,grp==g2 & keep_elec);
                        du_perm_post = coh_perm_POST(grp==g1 & keep_elec,grp==g2 & keep_elec);

                        dumm_mean_both{g1,g2}(:,fq) = nanmean([du ; du_post])';
                        dumm_perm_mean_both{g1,g2}(:,fq) = nanmean([du_perm ; du_perm_post])';
                    else
                        dumm_mean_both{g1,g2} = [];
                        dumm_perm_mean_both{g1,g2} = [];
                    end
                end
            end
        end

        if mk==1
            subdata{1,mk} = [subdata{1,mk} ; dumm_mean_both{3,1}];
            subdata_perm{1,mk} = [subdata_perm{1,mk} ; dumm_perm_mean_both{3,1}];
        else
            subdata{1,mk} = [subdata{1,mk} ; dumm_mean_both{2,1}];
            subdata_perm{1,mk} = [subdata_perm{1,mk} ; dumm_perm_mean_both{2,1}];
        end
    end
end

colors = [ 252 141 98 ; 102 194 165 ]/255;
colors_light = [ 252 189 164 ; 164 237 214 ]/255;

figure;set(gcf,'color',[1 1 1])
avg_FC = nanmean(subdata{1,1});
sem_FC = (nanstd(subdata{1,1})/sqrt(size(subdata{1,1},1)));
plot(conn.freq,avg_FC,'Color',colors(1,:),'LineWidth',2);hold on ;
ciplot(avg_FC-sem_FC,avg_FC+sem_FC,conn.freq,colors_light(1,:),0.5)
avg_FC = nanmean(subdata{1,2});
sem_FC = (nanstd(subdata{1,2})/sqrt(size(subdata{1,2},1)));
plot(conn.freq,avg_FC,'Color',colors(2,:),'LineWidth',2);hold on ;
ciplot(avg_FC-sem_FC,avg_FC+sem_FC,conn.freq,colors_light(2,:),0.5)
plot(conn.freq,nanmean(subdata_perm{1,1}),'--','Color',colors(1,:),'LineWidth',2);hold on ;
plot(conn.freq,nanmean(subdata_perm{1,2}),'--','Color',colors_light(2,:),'LineWidth',2)

xlim([0 30])
ylim([.02 .07])

pvals = []; zvals=[];
for fq = 1 : size(subdata{1,1},2)
    [pvals(1,fq),h,stats] = signrank(subdata{1,1}(:,fq),subdata_perm{1,1}(:,fq));
    zvals(1,fq) = stats.zval;
    [pvals(2,fq),h,stats] = signrank(subdata{1,2}(:,fq),subdata_perm{1,2}(:,fq));
    zvals(2,fq) = stats.zval;
end
[idx, idxs_mk1] = findenough(pvals(1,:),0.01,5,'<=');
[idx, idxs_mk2] = findenough(pvals(2,:),0.01,5,'<=');

plot(conn.freq(idxs_mk1),.023,'s','MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',6)
plot(conn.freq(idxs_mk2),.022,'s','MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'MarkerSize',6)

for mk = 1 : 2
    sub = nanmean(subdata{1,mk}(:,conn.freq>4 & conn.freq<50));
    subfq = conn.freq(conn.freq>4 & conn.freq<50);
    peakFQ(mk) = subfq(find(sub==max(sub)))
    minFQ = subfq(find(sub>=max(sub)-2*((max(sub)-min(sub))/3)));
    %minFQ = subfq(find(sub>=max(sub)-((max(sub)-min(sub))/2)));
    minmaxFQ(mk,:) = minFQ([1 end])

end
plot(peakFQ(1),.065,'.','Color',colors(1,:),'MarkerSize',30)
line([minmaxFQ(1,1) minmaxFQ(1,2)],[.065 .065],'Color',colors(1,:),'MarkerSize',30)
plot(peakFQ(2),.066,'.','Color',colors(2,:),'MarkerSize',30)
line([minmaxFQ(2,1) minmaxFQ(2,2)],[.066 .066],'Color',colors(2,:),'MarkerSize',30)
set(gca,"FontSize",16); xlabel('Frequency (Hz)');ylabel('Average coherency')
text(25,0.065,'mk L','Color', colors(1,:) ,"FontSize",16)
text(25,0.06,'mk H','Color', colors(2,:) ,"FontSize",16)

save([path2go 'rsSPK_foi.mat'],'peakFQ','minmaxFQ')

%% Run the true post processing using the extracted FQ of interest

clear
path2go = 'I:\Havoc_REC\';

contrast2do = {'Saline' 'DCZ'}

if strcmp(contrast2do{2},'DCZ')
    sessions = {'L082522a' 'L090122a' 'L081822b' 'L090822a' 'H041422a' 'H051222a' 'H042122a' 'H052422a' }
    cond = {'DCZ' 'DCZ'  'Saline' 'Saline' 'DCZ' 'DCZ' 'Saline' 'Saline' };
elseif strcmp(contrast2do{2},'CNO')
    sessions = {'L081822b' 'L090822a' 'L092222a' 'H042122a' 'H052422a' 'H060322a'}
    cond = {'Saline' 'Saline' 'CNO' 'Saline' 'Saline' 'CNO'};
end

load([path2go 'rsSPK_elecreject.mat']) %- load the rejection file
load([path2go 'rsSPK_foi.mat'],'peakFQ','minmaxFQ') %- load freq of interest

ref = 'bipol';
measure2take = 'coh'; %- can be 'powcorr', 'coh' or 'wpli_debiased'
data4glm=[];
for s = 1 : length(sessions)
    clearvars -except sessions cond ref path2go data4glm s measure2take peakFQ minmaxFQ elec_reject session_reject contrast2do
    list = dir([path2go sessions{s} '\*corr_permfilt_' ref '.mat']);
    if ~isempty(list)
        load([path2go sessions{s} '\' list(1).name])

        fq = conn.freq>=min(minmaxFQ(:,1)) & conn.freq<=max(minmaxFQ(:,2));

        eval(['coh = mean(conn(1).' measure2take '(:,:,fq),3);']);
        eval(['coh_POST = mean(conn_post.' measure2take '(:,:,fq),3);']);

        %- remove the diagonal (by adding a matrix of 0s with NaN on the diagonal!)
        v = NaN(length(coh),1);
        D = diag(v);
        coh = coh + D;
        coh_POST = coh_POST + D;
        if strcmp(sessions{s}(1),'L')
            grp = [ones(size(coh,1)/3,1) ; 2*ones(size(coh,1)/3,1) ; 3*ones(size(coh,1)/3,1)]; %- 2 elec with same contact number
            mk = 1;
        else
            grp = [ones(size(coh,1)/2,1) ; 2*ones(size(coh,1)/2,1)]; %- 2 elec with same contact number
            mk = 2 ;
        end
        %- remove elc that were identified as too noisy
        sess_rej = find(ismember(session_reject,sessions{s})); %- match current session with rejection matrix
        keep_elec = true(length(coh),1);
        keep_elec(elec_reject{sess_rej})=false;

        dumm_mean = [];
        dumm_std= [];
        dumm_post_mean = [];
        dumm_post_std= [];
        nbgrp = unique(grp);
        nbgrp(isnan(nbgrp))=[];

        for g1 = 1 : length(nbgrp)
            for g2 = 1 : length(nbgrp)
                if g1>=g2 %- only do that below diagonal, to remove the symetrical matrix
                    du = coh(grp==g1 & keep_elec,grp==g2 & keep_elec);
                    du_post = coh_POST(grp==g1 & keep_elec,grp==g2 & keep_elec);
                    dumm_mean(g1,g2) = nanmean(du(:));
                    dumm_std(g1,g2) = nanstd(du(:));
                    dumm_post_mean(g1,g2) = nanmean(du_post(:));
                    dumm_post_std(g1,g2) = nanstd(du_post(:));

                    data4glm = [data4glm ; nanmean(du_post-du)' , repmat(g1,size(du,2),1)  , repmat(g2,size(du,2),1)  , (1:size(du,2))' , repmat(mk,size(du,2),1) , repmat(s,size(du,2),1)];
                else
                    dumm_mean(g1,g2) = NaN;
                    dumm_std(g1,g2) = NaN;
                    dumm_post_mean(g1,g2) = NaN;
                    dumm_post_std(g1,g2) = NaN;
                end
            end
        end
    end

end

%- only consider the AMG-vlPFC for each monkey (ar3 + ar1 for Laf or 2 and 1 for Hav)
take = (data4glm(:,2)==3 & data4glm(:,3)==1 & data4glm(:,5)==1) | ...
    (data4glm(:,2)==2 & data4glm(:,3)==1 & data4glm(:,5)==2);

cond_names = cell(size(data4glm,1),1);
for i = 1: length(cond)
    cond_names(data4glm(:,end)==i,1)=cond(i);
end
data_tab = table(data4glm(take,1),cond_names(take),data4glm(take,4),data4glm(take,5),'VariableNames',{'coh','cond','elec','mk'});

data_tab.elec=categorical(data_tab.elec);
data_tab.mk=categorical(data_tab.mk);
data_tab.cond=categorical(data_tab.cond);

if strcmp(contrast2do{2},'DCZ')
    data_tab.cond = reordercats(data_tab.cond,{'Saline' 'DCZ' });
elseif strcmp(contrast2do{2},'CNO')
    data_tab.cond = reordercats(data_tab.cond,{'Saline' 'CNO' });
end


%- anova style
[p,tbl,stats,terms] = anovan(data_tab.coh,{data_tab.cond data_tab.mk data_tab.elec},'model',[1 0 0 ; 0 1 0 ; 1 1 0],...
                            'Varnames',{'cond' 'mk' 'elec'})
                        %  'random',3,'Varnames',{'cond' 'mk' 'elec'})
figure;multcompare(stats,'CType','bonferroni','Dimension',[1 2])
figure;multcompare(stats,'CType','bonferroni','Dimension',[1 ])

%- glmm style
lme_res = fitglme(data_tab, 'coh ~ 1 + cond*mk ')
lme_res = fitglme(data_tab, 'coh ~ 1 + cond*mk + (1|elec)','DummyVarCoding','effects')
anova(lme_res)

% emmeans(lme_res,{'cond' })

%- plot average coherency

colors = [  0.6510    0.8078    0.8902 ;
            0.1216    0.4706    0.7059 ;
            0.6980    0.8745    0.5412 ;
            0.2000    0.6275    0.1725 ;
            0.9843    0.6039    0.6000 ;
            0.8902    0.1020    0.1098]; %- originally using cbrewer but well, more stuff to have so skipped that!
if strcmp(contrast2do{2}(1:3),'CNO')
    colors = colors([1 2 5 6],:);
end

figure;
line([0 0],[0 length(contrast2do)+1],'Color','k');hold on
x = 1;
for cd = 1 : length(contrast2do)

    X= data_tab.coh(data_tab.cond==contrast2do{cd}) ;
    % X= (data_tab.coh(data_tab.mk==num2str(m) & data_tab.cond==cond{cd}) ;
    if ~isempty(X)
        yl=cd;
        wdth = .5;
        boxplot_ind(X,yl,wdth,[.4 .4 .4 ; colors(x+1,:)]);hold on
    end
    x = x+2;
end
set(gca,'view',[90 -90],'color','none','FontSize',16);
set(gca,'YTick',1:length(contrast2do),'YTickLabel',contrast2do,'YTickLabelRotation',25)
ylim([0 length(contrast2do)+1]);
xlim([-0.015 0.015])
xlabel('Coherency difference (POST - PRE)')

%% for revision: look at coherency in channels with increase/decrease of FR in vlPFC

clear
path2go = 'I:\Havoc_REC\';

contrast2do = {'DCZ' ''}

if strcmp(contrast2do{2},'DCZ')
    sessions = {'L082522a' 'L090122a' 'L081822b' 'L090822a' 'H041422a' 'H051222a' 'H042122a' 'H052422a' }
    cond = {'DCZ' 'DCZ'  'Saline' 'Saline' 'DCZ' 'DCZ' 'Saline' 'Saline' };
elseif strcmp(contrast2do{2},'CNO')
    sessions = {'L081822b' 'L090822a' 'L092222a' 'H042122a' 'H052422a' 'H060322a'}
    cond = {'Saline' 'Saline' 'CNO' 'Saline' 'Saline' 'CNO'};
else
    sessions = {'L082522a' 'L090122a' 'H041422a' 'H051222a' }
    cond = {'DCZ' 'DCZ'   'DCZ' 'DCZ' };
end

load([path2go 'rsSPK_elecreject.mat']) %- load the rejection file
load([path2go 'rsSPK_foi.mat'],'peakFQ','minmaxFQ') %- load freq of interest
load([path2go 'rsSPK_spk_inc_decr.mat']) %- load neurons increasing/decreasing FR

ref = 'bipol';
measure2take = 'coh'; %- can be 'powcorr', 'coh' or 'wpli_debiased'
data4glm=[];x=0;
for s = 1 : length(sessions)
    clearvars -except sessions cond ref path2go data4glm s measure2take peakFQ minmaxFQ elec_reject session_reject contrast2do spk_inc_decr change_coh x
    list = dir([path2go sessions{s} '\*corr_permfilt_' ref '.mat']);
    if ~isempty(list)
        load([path2go sessions{s} '\' list(1).name])

        fq = conn.freq>=min(minmaxFQ(:,1)) & conn.freq<=max(minmaxFQ(:,2));

        eval(['coh = mean(conn(1).' measure2take '(:,:,fq),3);']);
        eval(['coh_POST = mean(conn_post.' measure2take '(:,:,fq),3);']);

        %- remove the diagonal (by adding a matrix of 0s with NaN on the diagonal!)
        v = NaN(length(coh),1);
        D = diag(v);
        coh = coh + D;
        coh_POST = coh_POST + D;
        if strcmp(sessions{s}(1),'L')
            grp = [ones(size(coh,1)/3,1) ; 2*ones(size(coh,1)/3,1) ; 3*ones(size(coh,1)/3,1)]; %- 2 elec with same contact number
            mk = 1;
        else
            grp = [ones(size(coh,1)/2,1) ; 2*ones(size(coh,1)/2,1)]; %- 2 elec with same contact number
            mk = 2 ;
        end

        %- remove elc that were identified as too noisy
        sess_rej = find(ismember(session_reject,sessions{s})); %- match current session with rejection matrix
        keep_elec = true(length(coh),1);
        keep_elec(elec_reject{sess_rej})=false;

        nbgrp = unique(grp);
        nbgrp(isnan(nbgrp))=[];

        %- map the bipolar referencing
        if strcmp(sessions{s}(1),'H')
            bipolars = [(1:15)' , (2:16)' ;
                (17:31)' , (18:32)'];
        else
            bipolars = [(1:15)' , (2:16)' ;
                (17:31)' , (18:32)' ;
                (33:47)' , (34:48)'];
        end


        take_ch =[];
        dir_FR = [1 -1];
        for nn = 1 : 2
            for n = 1 : length(spk_inc_decr{nn})
                if ismember(sessions{s},{spk_inc_decr{nn}{n}(1:8)})~=0
                   take_ch = [take_ch ; dir_FR(nn) str2num(spk_inc_decr{nn}{n}(12:14))];
                end
            end
        end
        if ~isempty(take_ch)
            for g = 1 : length(take_ch(:,1))
                  [bip2consider,j] = find(bipolars==take_ch(g,2));
                  vlpfc_elec = zeros(size(keep_elec));
                  vlpfc_elec(bip2consider)=true;
                 % compare to grp1 elec (amygdala)
                    du = coh(grp==1 & keep_elec,vlpfc_elec & keep_elec);
                    du_post = coh_POST(grp==1 & keep_elec,vlpfc_elec & keep_elec);
 
                if ~isempty(du)
                     x = x + 1;
                    change_coh(x,:) = [nanmean(du_post(:)-du(:)) take_ch(g,1)];

                end
            end
        end
    end
end

%- Suppl Fig 4C
figure;
line([0 0],[0 length(contrast2do)+1],'Color','k');hold on
x = 1;
colors = [50 140 40 ; 50 200 40]/255
for cd = 1 : 2

    X= change_coh(change_coh(:,2)==dir_FR(cd),1) ;
    % X= (data_tab.coh(data_tab.mk==num2str(m) & data_tab.cond==cond{cd}) ;
    if ~isempty(X)
        yl=cd;
        wdth = .5;
        boxplot_ind(X,yl,wdth,[.4 .4 .4 ; colors(cd,:)]);hold on
    end
    x = x+2;
end
set(gca,'view',[90 -90],'color','none','FontSize',16);
set(gca,'YTick',1:length(contrast2do),'YTickLabel',{'vlPFC inc' 'vlPFC decr' },'YTickLabelRotation',25)
ylim([0 length(contrast2do)+1]);
xlim([-0.015 0.015])
xlabel('Coherency difference (POST - PRE)')

kruskalwallis(change_coh(:,1),change_coh(:,2))

