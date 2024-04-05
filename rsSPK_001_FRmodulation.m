%% Extract neurons' firing rate and compare PRE/POST injection
%- 
%- Need timestamps for every neurons + events!
%- Change line 87 "cond" to look at saline vs DCZ or Saline vs CNO
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2022.12

clear;

%- session to consider and matching condition
session = {'H041422a' 'H051222a' 'H042122a' 'H052422a' 'H060322a'...
    'L081822b' 'L082522a' 'L090122a' 'L090822a' 'L092222a'};
cond = { 'DCZ'    'DCZ' 'Saline' 'Saline' 'CNO' ...
    'Saline' 'DCZ' 'DCZ'    'Saline' 'CNO'};

pathname = 'I:\\Havoc_REC\'; %- file path
cd(pathname)

%- extract FR around events
all_unit = []; all_ch = []; all_evt = []; all_name={};
x = 0;
for ss = 1 : length(session)
    list_spk = dir([pathname 'SORTED NEURONS\' session{ss} '*Ch*Clus*'])
    sessionpath = [pathname session{ss} '\'];

    %- load event
    list_evt = dir([sessionpath '*KBD*.mat']);

    evts = [];
    for i = 1 : length(list_evt)
        load([sessionpath list_evt(i).name])
        if strcmp(session{ss},'H051222a') & i == 1% 2 mistake code in KBD2
            EVT = EVT(end);
        end
        evts = [evts , EVT];
    end
    evts = sortrows(evts'); %- sort them, normally 1 event in each, first one is the injection time, second is the start post-injection window

    %- process each spike
    for n = 1 : length(list_spk)
        x = x + 1;
        load([pathname 'SORTED NEURONS\' list_spk(n).name])

        %- extract psth around injection time (30sec bins)
        [psth_trials,~] = mpsth(timestamps,evts(1),'fr',1,'binsz',30*1000,'chart',0,'pre',15*60*1000,'post',10*60*1000);
        all_unit_pre(x,:) = psth_trials(:,2)';
        time_pre = psth_trials(:,1)/1000/60;

        %- extract psth around end injection time (30sec bins)
        [psth_trials,~] = mpsth(timestamps,evts(2),'fr',1,'binsz',30*1000,'chart',0,'pre',10*60*1000,'post',20*60*1000);
        time_post = psth_trials(:,1)/1000/60;
        all_unit_post(x,:) = psth_trials(:,2)';

        %- time vector with the multiple alignments
        time = [time_pre' ,time_post' ; ones(size(time_pre))' , 2*ones(size(time_post))'];

        all_ch(x) = str2num(list_spk(n).name(12:14)); %- elec contact
        all_evt(x,:) = evts; %- event times
        all_cond{x} = cond(ss); %- injection type
        all_mk{x} = list_spk(n).name(1); %- monkey name
        all_name{x} = list_spk(n).name; %- monkey name

        %- extract area info (based on monkey name and elec contact)
        if strcmp(all_mk{x},'H') && all_ch(x)<=16
            all_area{x} = 'AMG';
        elseif strcmp(all_mk{x},'H') && all_ch(x)>16
            all_area{x} = 'vlPFC';
        elseif strcmp(all_mk{x},'L') && all_ch(x)<=32
            all_area{x} = 'AMG';
        elseif strcmp(all_mk{x},'L') && all_ch(x)>32
            all_area{x} = 'vlPFC';
        end

    end

end
all_unit = [all_unit_pre all_unit_post];

%- normalization (referenced to the pre-injection time)
clear all_norm
ref = time(1,:)>=-10 & time(1,:)<0 & time(2,:)==1;
all_norm = (all_unit - repmat(mean(all_unit(:,ref),2),1,size(all_unit,2))) ./ repmat(std(all_unit(:,ref)')',1,size(all_unit,2));

%% Figure + data extraction : Average normalized FR across condition and areas

cond = {'Saline' 'DCZ'} %- condition to look at
arr = {'AMG' 'vlPFC'} %- area to look at
nb_bins_thr = 4; %- number of bins above ~2sd to be considered sig

colors = [  0.6510    0.8078    0.8902 ;
            0.1216    0.4706    0.7059 ;
            0.6980    0.8745    0.5412 ;
            0.2000    0.6275    0.1725 ;
            0.9843    0.6039    0.6000 ;
            0.8902    0.1020    0.1098]; %- originally using cbrewer but well, more stuff to have so skipped that! 
if strcmp(cond{2}(1:3),'CNO')
    colors = colors([1 2 5 6],:);
end

%- plot Fig 5B,C
figure;
clear prepost
x1 = []; y_ar = []; y_mk = []; y_cd = [];
for ar = 1 : length(arr)
    x=1;
    subplot(length(arr),1,ar);
    for cd = 1 : length(cond)

        %- take only the current cond + area
        takeme = ismember(all_area,arr{ar}) & ismember([all_cond{:}],cond{cd})    ;
        takeme_mk1 = ismember(all_area,arr{ar}) & ismember([all_cond{:}],cond{cd}) & ismember(all_mk,'L')   ; %- same for each monkey
        takeme_mk2 = ismember(all_area,arr{ar}) & ismember([all_cond{:}],cond{cd}) & ismember(all_mk,'H')   ;

        %   plot(  mean(all_norm(takeme,:)),'Color',colors(x,:)); hold on

        avg_FR = smooth( mean(all_norm(takeme,:)),4,'moving'); %- smoothing across 4 bins = 2min
        sem_FR = smooth( std(all_norm(takeme,:))./sqrt(sum(takeme)),4,'moving');

       % plot( avg_FR( (time(1,:)<0 & time(2,:)==1)|(time(1,:)>0 & time(2,:)==2))  ,'Color',colors(x+1,:),'LineWidth',2); hold on
        plot( avg_FR,'Color',colors(x+1,:),'LineWidth',2); hold on
        ciplot(avg_FR-sem_FR,avg_FR+sem_FR,1:length(avg_FR),colors(x,:),.4)

        % plot(  smooth(  mean(all_norm(takeme,:)) ,4,'moving')  ,'Color',colors(x+1,:),'LineWidth',2);hold on
        % plot(  smooth(  mean(all_norm(takeme_mk1,:)) ,4,'moving')  ,'Color',colors(x,:),'LineWidth',1)
        % plot(  smooth(  mean(all_norm(takeme_mk2,:)) ,4,'moving')  ,'Color',colors(x,:),'LineWidth',1)

        text(3,0.5-(x/9),[cond{cd} ' - n=' num2str(sum(takeme))],'Color',colors(x+1,:),'FontSize',16)
        x=x+2;

        %- extract average zscore for PRE/POST across conditions
        prepost{ar,cd} = [mean(all_norm(takeme,time(1,:)<0 & time(2,:)==1),2) , ...
            mean(all_norm(takeme,time(1,:)>0 & time(2,:)==2),2)];

        %- find neurons with strong change in FR
        clear sig
        dumm = abs(all_norm(takeme,time(1,:)>0 & time(2,:)==2))>1.96;
        %         for i = 1 : size(dumm,1)
        %             [idx,idxs] = findenough(dumm(i,:),1,2,'==');
        %             if ~isempty(idx)
        %                 sig(i) = 1;
        %             else
        %                 sig(i)=0;
        %             end
        %         end

        sig = (sum(dumm,2)>=nb_bins_thr)';

        %- count it
        nSig_bin(ar,cd)=sum(sig);
        nTot_bin(ar,cd)=length(sig);

        %- put that into a matrix needed for binomial glm
        nUnits = size(sig',1);
        x1 = [x1 ; logical(sig)' ];
        y_ar = [y_ar ; repmat(arr{ar}(1:3),nUnits,1)];
        y_mk = [y_mk ; all_mk(takeme)'];
        y_cd = [y_cd ; repmat(cond{cd}(1:3),nUnits,1)];
    end

    %- some leftover figure stuff
    lim = find(time(2,:)==2,1,'first')+0.5;
    t0 =  find(time(1,:)>=0,1,'first'); t1 =  find(time(1,:)<=0,1,'last');
    line([t0 t0],[-1.5 1.5],'Color','k','LineWidth',1)
    line([t1 t1],[-1.5 1.5],'Color','k','LineWidth',1)
    line([0 length(all_norm(1,:))+1],[0 0],'Color','k','LineWidth',1)
    line([lim lim],[-2.5 2.5],'Color','w','LineWidth',3)
    set(gca,'Xtick',1:10:length(time),'Xticklabel',time(1,1:10:end),'FontSize',16)
    xlim([0 length(time(1,time(1,:)<=20))+1])
    ylim([-0.5 1])
    title(arr{ar})
    xlabel('Time (min)')
    ylabel('Average normalized firing rate (sd)')
end


%- binomial model for sig neurons
modeldata = table(x1,y_ar,y_mk,y_cd,'VariableNames',{'sig' 'area' 'mk' 'cd'});
lme_percSig = fitglme(modeldata,'sig ~ 1 + area*cd + (1|mk) ','Distribution','binomial','Link','logit');

anova(lme_percSig)

res = grpstats(modeldata,{'area' 'mk' 'cd'},"mean")

%- plot % modulated neurons (e.g. Fig 6E)
figure;
ax = bar(100*(nSig_bin./nTot_bin));hold on
ax(1).FaceColor = colors(2,:); ax(1).DisplayName = cond{1};
ax(2).FaceColor = colors(4,:); ax(2).DisplayName = cond{2};
markers4mk = 'ov';
xax = [.85 1.15 1.85 2.15 ];
mkk={'L' 'H'};
for mk = 1 : 2
    plot(xax,100*res(ismember(res.mk,mkk{mk}),:).mean_sig,markers4mk(mk),'MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k','DisplayName',['mk ' mkk{mk}])
end
set(gca,'XTick',1:length(arr),'XTickLabel',arr,'FontSize',16)
ylabel('Proportion of modulated neurons')
legend("show",'Location','northeastoutside')


%% Plot for each neuron the difference in PRE/POST FR (Fig 5B,C)
figure;
for ar = 1 : length(arr)
    x = 1;
    subplot(1,length(arr),ar) ;line([0 0],[0 length(cond)+1],'Color','k');title(arr{ar});hold on
    for cd = 1 : length(cond)

        X= prepost{ar,cd}(:,2)-prepost{ar,cd}(:,1) ;
        if ~isempty(X)
            yl=cd;
            wdth = .5;

            boxplot_ind(X,yl,wdth,[.4 .4 .4 ; colors(x+1,:)]);hold on
            [p_wilcox(ar,cd),~,sta] = signrank(prepost{ar,cd}(:,1),prepost{ar,cd}(:,2));
            text(max(X)+1.5,cd,['p=' num2str(round(p_wilcox(ar,cd)*1000)/1000)],'HorizontalAlignment','center')
        end
        x = x+2;

    end
    set(gca,'view',[90 -90],'color','none','FontSize',16);
    set(gca,'YTick',1:length(cond),'YTickLabel',cond,'YTickLabelRotation',25)
    ylim([0 length(cond)+1]);
    xlim([-10 10])
    xlabel('Normalized firing rate difference (POST - PRE) (sd)')
end

%- compare the distributions between the 2 conditions
for ar = 1 : length(arr)
    X = prepost{ar,1}(:,2)-prepost{ar,1}(:,1) ;
    Y = prepost{ar,2}(:,2)-prepost{ar,2}(:,1) ;
    data = [X , ones(size(X)) ; Y , 2*ones(size(Y))];
    [p_kw_vehdcz(ar),tbl_kw_vehdcz{ar},sta] = kruskalwallis(data(:,1),data(:,2),'off');
end


%% for revision: extract neurons in vlPFC that increase or decrease FR
%- this will be used to split channels and look for difference in coherency
%- depending on whether neurons on these channels increase/decrease FR

%- take only the current cond + area
takeme = ismember(all_area,'vlPFC') & ismember([all_cond{:}],'DCZ')    ;

dumm = [mean(all_norm(takeme,time(1,:)<0 & time(2,:)==1),2) , ...
        mean(all_norm(takeme,time(1,:)>0 & time(2,:)==2),2)];

name_sub = all_name(takeme)';

spk_inc_decr{1} = name_sub(dumm(:,2)-dumm(:,1)>0);
spk_inc_decr{2} = name_sub(dumm(:,2)-dumm(:,1)<0);

save('rsSPK_spk_inc_decr.mat','spk_inc_decr')

%% imagesc Fig S4A,B

figure;       
clear prepost
for ar = 1 : length(arr)
    x=1;
    subplot(1,length(arr),ar);
    all_plot =[];lim_pl =[];
    for cd = 2 %: length(cond)
        takeme = ismember(all_area,arr{ar}) & ismember([all_cond{:}],cond{cd})    ;

        dumm = all_norm(takeme,:);
        avg_end = mean(dumm(:,end-20:end),2);
        [~,d] = sortrows(avg_end);

        all_plot = [all_plot ; dumm(d,:)];
        lim_pl = [lim_pl ; length(all_plot(:,1))+.5];

    end

    imagesc(all_plot,[-1.5 1.5]);ylabel('neurons')
    hold on
    for j = 1 : length(lim_pl)-1
        line([0 length(all_plot)],[lim_pl(j) lim_pl(j)],'Color','k')
    end
    set(gca,'Xtick',1:10:length(time),'Xticklabel',time(1,1:10:end),'FontSize',16)
    lim = find(time(2,:)==2,1,'first')+0.5;
    line([lim lim],[0 length(all_plot)],'Color','w','LineWidth',3)
    t0 =  find(time(1,:)>=0,1,'first'); t1 =  find(time(1,:)<=0,1,'last');
    line([t0 t0],[0 length(all_plot)],'Color','k','LineWidth',1)
    line([t1 t1],[0 length(all_plot)],'Color','k','LineWidth',1)

    title(arr{ar})

end
   
%% to give data to Cat:
%    anova_tab =[];
%    for ar = 1  : size(prepost,1)
%        for cd = 1 : 3
%              anova_tab = [anova_tab ; prepost{ar,cd}(:,1) prepost{ar,cd}(:,2) ar*ones(size( prepost{ar,cd}(:,1))) cd*ones(size( prepost{ar,cd}(:,1)))];
%        end
%    end
%    anova_tab_header = {'PRE norm FR' 'POST norm FR' 'area' 'cond'}
%    save('I:\Havoc_REC\rsSPK_normFR.mat','anova_tab','anova_tab_header')
   
