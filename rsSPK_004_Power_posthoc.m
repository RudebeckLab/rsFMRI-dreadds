%% LFP Power post-processing
%-
%- Require fieldtrip and to have run "..._002_LFPpreproc.m"
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.01

clear

path2go = 'I:\Havoc_REC\'; %- location of the folders with the data
conds = {'Saline' 'CNO'}; %- takes only 2 inputs

%- check how much channel rejected
load([path2go 'rsSPK_elecreject.mat']) %- load the rejection file

nb_rej=[];tot_ch=[];
for i = 1 : length(elec_reject)
    nb_rej(i)=size(elec_reject{i},2);
    if strcmp(session_reject{i}(1),'L')
        tot_ch(i)=45;
    else
        tot_ch(i)=30;
    end
end
mean(nb_rej./tot_ch)

%- take the sessions according to your choice
if strcmp(conds{2},'DCZ')
    sessions = {'L082522a' 'L090122a' 'L081822b' 'L090822a' 'H041422a' 'H051222a' 'H042122a' 'H052422a' }
    cond = {'DCZ' 'DCZ'  'Saline' 'Saline' 'DCZ' 'DCZ' 'Saline' 'Saline' };
elseif strcmp(conds{2},'CNO')
    sessions = {'L081822b' 'L090822a' 'L092222a' 'H042122a' 'H052422a' 'H060322a'}
    cond = {'Saline' 'Saline' 'CNO' 'Saline' 'Saline' 'CNO'};
end

%- which ref / measure
ref = 'bipol';
measure2take = {'coh' 'coh_perm'}; %- can be 'powcorr', 'coh' or 'wpli_debiased'

figure;set(gcf,'color',[1 1 1])

all_elec = [];
all_mk = [];
all_cond = [];
all_fft = [];
all_fft_pre_raw =[];
all_fft_post_raw =[];
for s = 1 : length(sessions)
    list = dir([path2go sessions{s} '\*corr_permfilt_' ref '.mat']);
    if ~isempty(list)
        load([path2go sessions{s} '\' list(1).name])

        %- define monkey and elec groups
        if strcmp(sessions{s}(1),'L')
            grp = [ones(size(conn.label,1)/3,1) ; 2*ones(size(conn.label,1)/3,1) ; 3*ones(size(conn.label,1)/3,1)]; %- 2 elec with same contact number
            mk = 1;
        else
            grp = [ones(size(conn.label,1)/2,1) ; 2*ones(size(conn.label,1)/2,1)];
            mk = 2 ;
        end

        %- remove elc that were identified as too noisy
        sess_rej = find(ismember(session_reject,sessions{s})); %- match current session with rejection matrix
        keep_elec = true(length(conn.label),1);
        keep_elec(elec_reject{sess_rej})=false;

        % extract elec groups
        nbgrp = unique(grp);
        nbgrp(isnan(nbgrp))=[];

        %- compute normalized power for post
        fft_post.powspctrm_norm = (fft_post.powspctrm - fft_pre.powspctrm) ./ fft_pre.powspctrm_std ;

        %- extract all useful info
        all_fft = [all_fft ; fft_post.powspctrm_norm(keep_elec,:)];
        all_fft_pre_raw = [all_fft_pre_raw ; fft_pre.powspctrm(keep_elec,:)];
        all_fft_post_raw = [all_fft_post_raw ; fft_post.powspctrm(keep_elec,:)];
        all_elec = [all_elec ; grp(keep_elec) ];
        all_mk = [all_mk ; repmat(mk,size(grp(keep_elec))) ];
        all_cond = [all_cond ; repmat(cond(s),size(grp(keep_elec))) ];

        %- some averages
        fft_all =[];
        fft_sem =[];
        for g1 = 1 : length(nbgrp)
            fft_all(g1,:) = nanmean(fft_post.powspctrm_norm(grp==g1 & keep_elec,:));
            fft_sem(g1,:) = nanstd(fft_post.powspctrm_norm(grp==g1 & keep_elec,:))/sqrt(sum(grp==g1 & keep_elec));
        end

        %-quick and dirty plot, per session
        subplot(2,4,s);
        ciplot(fft_all-fft_sem,fft_all+fft_sem,conn.freq);hold on
        plot(conn.freq,fft_all');hold on
        xlim([1 30])

    end
end

%- colors for plot
colors = [ 0.1216    0.4706    0.7059 ;
    0.2000    0.6275    0.1725 ;
    0.8902    0.1020    0.1098]; %- originally using cbrewer but well, more stuff to have in the folder so skipped that!
if strcmp(conds{2}(1:3),'CNO')
    colors = colors([1 3],:);
end

figure;
clear takeme_amg takeme_vlpfc

for cd = 1 : length(conds)
    takeme_amg{cd} =  (all_mk==1 & all_elec==1 & ismember(all_cond,conds{cd})) | (all_mk==2 & all_elec==1 & ismember(all_cond,conds{cd}));
    takeme_vlpfc{cd} =  (all_mk==1 & all_elec==3 & ismember(all_cond,conds{cd})) | (all_mk==2 & all_elec==2 & ismember(all_cond,conds{cd}));

    avg_pow = nanmean(all_fft(takeme_amg{cd},:));
    sem_pow = (nanstd(all_fft(takeme_amg{cd},:))/sqrt(sum(takeme_amg{cd})));
    subplot(1,2,1);
    line([0 100],[0 0],'Color','k')
    plot(conn.freq,avg_pow,'Color',colors(cd,:),'LineWidth',2);hold on ;
    ciplot(avg_pow-sem_pow,avg_pow+sem_pow,conn.freq,colors(cd,:),0.5)

    pvals = []; zvals=[];
    for fq = 1 : size(conn.freq,2)
        [pvals(1,fq),h,stats] = signrank(all_fft(takeme_amg{cd},fq));
        zvals(1,fq) = stats.zval;
    end
    [idx, idxs] = findenough(pvals(1,:),0.01,5,'<=');
    if ~isempty(idxs)
        plot(conn.freq(idxs),-2.5+(cd)/10,'s','MarkerEdgeColor',colors(cd,:),'MarkerFaceColor',colors(cd,:),'MarkerSize',6)
    end

    xlim([0 100]);ylim([-3 3])
    set(gca,"FontSize",16) ; xlabel('Frequency (Hz)') ; ylabel('Normalized Power (sd)'),

    if cd == 2
        pvals = []; zvals=[];
        for fq = 1 : size(conn.freq,2)
            dumm = [all_fft(takeme_amg{1},fq) , ones(size(all_fft(takeme_amg{1},fq))) ;
                all_fft(takeme_amg{2},fq) , 2*ones(size(all_fft(takeme_amg{2},fq)))];
            [pvals(1,fq),h,stats] = kruskalwallis(dumm(:,1),dumm(:,2),'off');
        end
        [idx, idxs] = findenough(pvals(1,:),0.01,5,'<=');
        if ~isempty(idxs)

            plot(conn.freq(idxs),-2.5+(cd+1)/10,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
        end
    end

    %- same for vlpfc
    subplot(1,2,2);
    avg_pow = nanmean(all_fft(takeme_vlpfc{cd},:));
    sem_pow = (nanstd(all_fft(takeme_vlpfc{cd},:))/sqrt(sum(takeme_vlpfc{cd})));
    line([0 100],[0 0],'Color','k')
    plot(conn.freq,avg_pow,'Color',colors(cd,:),'LineWidth',2);hold on ;
    ciplot(avg_pow-sem_pow,avg_pow+sem_pow,conn.freq,colors(cd,:),0.5)

    pvals = []; zvals=[];
    for fq = 1 :  size(conn.freq,2)
        [pvals(1,fq),h,stats] = signrank(all_fft(takeme_vlpfc{cd},fq));
        zvals(1,fq) = stats.zval;
    end
    [idx, idxs] = findenough(pvals(1,:),0.01,5,'<=');
    if ~isempty(idxs)
        plot(conn.freq(idxs),-2.5+cd/10,'s','MarkerEdgeColor',colors(cd,:),'MarkerFaceColor',colors(cd,:),'MarkerSize',6)
    end
    xlim([0 100]);ylim([-3 3])
    set(gca,"FontSize",16) ; xlabel('Frequency (Hz)') ; ylabel('Normalized Power (sd)'),

    if cd == 2
        pvals = []; zvals=[];
        for fq = 1 :  size(conn.freq,2)
            dumm = [all_fft(takeme_vlpfc{1},fq) , ones(size(all_fft(takeme_vlpfc{1},fq))) ;
                all_fft(takeme_vlpfc{2},fq) , 2*ones(size(all_fft(takeme_vlpfc{2},fq)))];
            [pvals(1,fq),h,stats] = kruskalwallis(dumm(:,1),dumm(:,2),'off');
        end
        [idx, idxs] = findenough(pvals(1,:),0.01,5,'<=');
        if ~isempty(idxs)
            plot(conn.freq(idxs),-2.5+(cd+1)/10,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
        end
    end

end
