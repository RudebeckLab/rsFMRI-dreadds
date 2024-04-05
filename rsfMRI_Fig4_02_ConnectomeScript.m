%Load in data
clear all; close all;
cd C:/Users/Catherine/Documents/Data/rsMRI/
CHARMData = load('CHARMData.mat');
CHARMData = CHARMData.ans;
SARMData = load('SARMData.mat');
SARMData = SARMData.ans;
D99Data = load('D99Data.mat');
D99Data = D99Data.ans;
CHARM_MAIN = CHARMData;
SARM_MAIN = SARMData;
D99_MAIN = D99Data;

%% Prepare the data
%split data
charm_zval = CHARM_MAIN{2,1};
charm_sbj = CHARM_MAIN{2,2};
charm_drug = CHARM_MAIN{2,3};
charm_sess = CHARM_MAIN{2,4};
charm_connec = CHARM_MAIN{2,5};
sarm_zval = SARM_MAIN{2,1};
sarm_sbj = SARM_MAIN{2,2};
sarm_drug = SARM_MAIN{2,3};
sarm_sess = SARM_MAIN{2,4};
sarm_connec = SARM_MAIN{2,5};
d99_zval = D99_MAIN{2,1};
d99_sbj = D99_MAIN{2,2};
d99_drug = D99_MAIN{2,3};
d99_sess = D99_MAIN{2,4};
d99_connec = D99_MAIN{2,5};

drug_num = unique(charm_drug);
sess_num = unique(charm_sess);
sbj_num = unique(charm_sbj);

%analyze average FC for each subject
charm_drugdata = []; sarm_drugdata = []; d99_drugdata = [];
for i = 1:size(drug_num,1)
    for k = 1:size(sbj_num,1)
        charm_drugdata(i,k) = mean(charm_zval(charm_drug==i & charm_sbj==k));
        sarm_drugdata(i,k) = mean(sarm_zval(sarm_drug==i & sarm_sbj==k));
        d99_drugdata(i,k) = mean(d99_zval(d99_drug==i & d99_sbj==k));
    end
end

for i=1:3
    for l=1:2
        charm_by_drug_and_subj_std(i,l) = std(charm_zval(charm_drug==i & charm_sbj==l),'omitnan');
        sarm_by_drug_and_subj_std(i,l) = std(sarm_zval(sarm_drug==i & sarm_sbj==l),'omitnan');
        d99_by_drug_and_subj_std(i,l) = std(d99_zval(d99_drug==i & d99_sbj==l),'omitnan');
    end
end


for i=1:3
    for l=1:2
        charm_by_drug_and_subj_sem(i,l) = charm_by_drug_and_subj_std(i,l)/(length(charm_zval(charm_drug==i & charm_sbj==l)).^0.5);
        sarm_by_drug_and_subj_sem(i,l) = sarm_by_drug_and_subj_std(i,l)/(length(sarm_zval(sarm_drug==i & sarm_sbj==l)).^0.5);
        d99_by_drug_and_subj_sem(i,l) = d99_by_drug_and_subj_std(i,l)/(length(d99_zval(d99_drug==i & d99_sbj==l)).^0.5);
    end
end

%Separate out data from each session
charm_veh1_hav = charm_zval(charm_drug==1 & charm_sess==1 & charm_sbj==1);
charm_veh2_hav = charm_zval(charm_drug==1 & charm_sess==2 & charm_sbj==1);
charm_dcz1_hav = charm_zval(charm_drug==2 & charm_sess==1 & charm_sbj==1);
charm_dcz2_hav = charm_zval(charm_drug==2 & charm_sess==2 & charm_sbj==1);
charm_cno1_hav = charm_zval(charm_drug==3 & charm_sess==1 & charm_sbj==1);
charm_veh1_laf = charm_zval(charm_drug==1 & charm_sess==1 & charm_sbj==2);
charm_veh2_laf = charm_zval(charm_drug==1 & charm_sess==2 & charm_sbj==2);
charm_dcz1_laf = charm_zval(charm_drug==2 & charm_sess==1 & charm_sbj==2);
charm_dcz2_laf = charm_zval(charm_drug==2 & charm_sess==2 & charm_sbj==2);
charm_cno1_laf = charm_zval(charm_drug==3 & charm_sess==1 & charm_sbj==2);

charm_veh_hav = ((charm_veh1_hav+charm_veh2_hav)/2);
charm_dcz_hav = ((charm_dcz1_hav+charm_dcz2_hav)/2);
charm_cno_hav = charm_cno1_hav;
charm_veh_laf = ((charm_veh1_laf+charm_veh2_laf)/2);
charm_dcz_laf = ((charm_dcz1_laf+charm_dcz2_laf)/2);
charm_cno_laf = charm_cno1_laf;

sarm_veh1_hav = sarm_zval(sarm_drug==1 & sarm_sess==1 & sarm_sbj==1);
sarm_veh2_hav = sarm_zval(sarm_drug==1 & sarm_sess==2 & sarm_sbj==1);
sarm_dcz1_hav = sarm_zval(sarm_drug==2 & sarm_sess==1 & sarm_sbj==1);
sarm_dcz2_hav = sarm_zval(sarm_drug==2 & sarm_sess==2 & sarm_sbj==1);
sarm_cno1_hav = sarm_zval(sarm_drug==3 & sarm_sess==1 & sarm_sbj==1);
sarm_veh1_laf = sarm_zval(sarm_drug==1 & sarm_sess==1 & sarm_sbj==2);
sarm_veh2_laf = sarm_zval(sarm_drug==1 & sarm_sess==2 & sarm_sbj==2);
sarm_dcz1_laf = sarm_zval(sarm_drug==2 & sarm_sess==1 & sarm_sbj==2);
sarm_dcz2_laf = sarm_zval(sarm_drug==2 & sarm_sess==2 & sarm_sbj==2);
sarm_cno1_laf = sarm_zval(sarm_drug==3 & sarm_sess==1 & sarm_sbj==2);

sarm_veh_hav = ((sarm_veh1_hav+sarm_veh2_hav)/2);
sarm_dcz_hav = ((sarm_dcz1_hav+sarm_dcz2_hav)/2);
sarm_cno_hav = sarm_cno1_hav;
sarm_veh_laf = ((sarm_veh1_laf+sarm_veh2_laf)/2);
sarm_dcz_laf = ((sarm_dcz1_laf+sarm_dcz2_laf)/2);
sarm_cno_laf = sarm_cno1_laf;

d99_veh1_hav = d99_zval(d99_drug==1 & d99_sess==1 & d99_sbj==1);
d99_veh2_hav = d99_zval(d99_drug==1 & d99_sess==2 & d99_sbj==1);
d99_dcz1_hav = d99_zval(d99_drug==2 & d99_sess==1 & d99_sbj==1);
d99_dcz2_hav = d99_zval(d99_drug==2 & d99_sess==2 & d99_sbj==1);
d99_cno1_hav = d99_zval(d99_drug==3 & d99_sess==1 & d99_sbj==1);
d99_veh1_laf = d99_zval(d99_drug==1 & d99_sess==1 & d99_sbj==2);
d99_veh2_laf = d99_zval(d99_drug==1 & d99_sess==2 & d99_sbj==2);
d99_dcz1_laf = d99_zval(d99_drug==2 & d99_sess==1 & d99_sbj==2);
d99_dcz2_laf = d99_zval(d99_drug==2 & d99_sess==2 & d99_sbj==2);
d99_cno1_laf = d99_zval(d99_drug==3 & d99_sess==1 & d99_sbj==2);

d99_veh_hav = ((d99_veh1_hav+d99_veh2_hav)/2);
d99_dcz_hav = ((d99_dcz1_hav+d99_dcz2_hav)/2);
d99_cno_hav = d99_cno1_hav;
d99_veh_laf = ((d99_veh1_laf+d99_veh2_laf)/2);
d99_dcz_laf = ((d99_dcz1_laf+d99_dcz2_laf)/2);
d99_cno_laf = d99_cno1_laf;

charm_sessdata = []; sarm_sessdata = []; d99_sessdata = [];
%Put the connectome data in vector format for ANOVA
charm_sessdata = [mean(charm_veh1_hav) mean(charm_veh1_laf);
    mean(charm_veh2_hav) mean(charm_veh2_laf); 
    mean(charm_dcz1_hav) mean(charm_dcz1_laf); 
    mean(charm_dcz2_hav) mean(charm_dcz2_laf);
    mean(charm_cno1_hav) mean(charm_cno1_laf)];

sarm_sessdata = [mean(sarm_veh1_hav) mean(sarm_veh1_laf);
    mean(sarm_veh2_hav) mean(sarm_veh2_laf); 
    mean(sarm_dcz1_hav) mean(sarm_dcz1_laf); 
    mean(sarm_dcz2_hav) mean(sarm_dcz2_laf);
    mean(sarm_cno1_hav) mean(sarm_cno1_laf)];

d99_sessdata = [mean(d99_veh1_hav) mean(d99_veh1_laf);
    mean(d99_veh2_hav) mean(d99_veh2_laf); 
    mean(d99_dcz1_hav) mean(d99_dcz1_laf); 
    mean(d99_dcz2_hav) mean(d99_dcz2_laf);
    mean(d99_cno1_hav) mean(d99_cno1_laf)];
    
charm_by_drug_mean = nanmean(charm_drugdata,2);
sarm_by_drug_mean = nanmean(sarm_drugdata,2);
d99_by_drug_mean = nanmean(d99_drugdata,2);
for i=1:3
    charm_by_drug_std(i) = nanstd(charm_drugdata(i,:));
    sarm_by_drug_std(i) = nanstd(sarm_drugdata(i,:));
    d99_by_drug_std(i) = nanstd(d99_drugdata(i,:));
end
for i=1:3
    charm_by_drug_sem(i) = charm_by_drug_std(:,i)/(sum(~isnan(charm_drugdata(i,:)),2).^0.5);
    sarm_by_drug_sem(i) = sarm_by_drug_std(:,i)/(sum(~isnan(sarm_drugdata(i,:)),2).^0.5);
    d99_by_drug_sem(i) = d99_by_drug_std(:,i)/(sum(~isnan(d99_drugdata(i,:)),2).^0.5);
end

charm_by_sess_mean = nanmean(charm_sessdata,2);
sarm_by_sess_mean = nanmean(sarm_sessdata,2);
d99_by_sess_mean = nanmean(d99_sessdata,2);
for i=1:5
    charm_by_sess_sem(i) = nanstd(charm_sessdata(i,:)')/(sum(~isnan(charm_sessdata(i,:)),2).^0.5);
    sarm_by_sess_sem(i) = nanstd(sarm_sessdata(i,:)')/(sum(~isnan(sarm_sessdata(i,:)),2).^0.5);
    d99_by_sess_sem(i) = nanstd(d99_sessdata(i,:)')/(sum(~isnan(d99_sessdata(i,:)),2).^0.5);
end

charmlength = length(charm_veh_hav);
sarmlength = length(sarm_veh_hav);
d99length = length(d99_veh_hav);

sesscharmveh(1:charmlength,1) = 1;
sesscharmdcz(1:charmlength,1) = 2;
sesscharmcno(1:charmlength,1) = 3;
sesscharmses1(1:charmlength,1) = 1;
sesscharmses2(1:charmlength,1) = 2;
sesscharmhav(1:charmlength,1) = 1;
sesscharmlaf(1:charmlength,1) = 2;
sesscharmconnec = charm_connec(1:charmlength,1);

sesssarmveh(1:sarmlength,1) = 1;
sesssarmdcz(1:sarmlength,1) = 2;
sesssarmcno(1:sarmlength,1) = 3;
sesssarmses1(1:sarmlength,1) = 1;
sesssarmses2(1:sarmlength,1) = 2;
sesssarmhav(1:sarmlength,1) = 1;
sesssarmlaf(1:sarmlength,1) = 2;
sesssarmconnec = sarm_connec(1:561,1);

sessd99veh(1:d99length,1) = 1;
sessd99dcz(1:d99length,1) = 2;
sessd99cno(1:d99length,1) = 3;
sessd99ses1(1:d99length,1) = 1;
sessd99ses2(1:d99length,1) = 2;
sessd99hav(1:d99length,1) = 1;
sessd99laf(1:d99length,1) = 2;
sessd99connec = d99_connec(1:d99length,1);


charm_sessgroup = [charm_sessdata(1,1) charm_sessdata(2,1);
    charm_sessdata(3,1) charm_sessdata(4,1);
    charm_sessdata(5,1) NaN;
    charm_sessdata(1,2) charm_sessdata(2,2);
    charm_sessdata(3,2) charm_sessdata(4,2);
    charm_sessdata(5,2) NaN];

sarm_sessgroup = [sarm_sessdata(1,1) sarm_sessdata(2,1);
    sarm_sessdata(3,1) sarm_sessdata(4,1);
    sarm_sessdata(5,1) NaN;
    sarm_sessdata(1,2) sarm_sessdata(2,2);
    sarm_sessdata(3,2) sarm_sessdata(4,2);
    sarm_sessdata(5,2) NaN];

d99_sessgroup = [d99_sessdata(1,1) d99_sessdata(2,1);
    d99_sessdata(3,1) d99_sessdata(4,1);
    d99_sessdata(5,1) NaN;
    d99_sessdata(1,2) d99_sessdata(2,2);
    d99_sessdata(3,2) d99_sessdata(4,2);
    d99_sessdata(5,2) NaN];

%% Compare VEH and DCZ treatments
%Let's organize our ANOVAn as: Drug, Session, Subject, Connection. So level one will
%be 1= VEH (H then L), 2 = DCZ. [1 1 1 1 2 2 2 2]. 
%Session will be in order 1-2 [1 2 1 2 1 2 1 2]. 
%Subject is H = 1, L = 2. [1 1 2 2 1 1 2 2].
%Connnection matrices are made with the x_connec variables.


sessalld99_nocno = [d99_veh1_hav;d99_veh2_hav;d99_veh1_laf;d99_veh2_laf; ...
    d99_dcz1_hav;d99_dcz2_hav;d99_dcz1_laf;d99_dcz2_laf];
sessallcharm_nocno = [charm_veh1_hav;charm_veh2_hav;charm_veh1_laf;charm_veh2_laf; ...
    charm_dcz1_hav;charm_dcz2_hav;charm_dcz1_laf;charm_dcz2_laf];
sessallsarm_nocno = [sarm_veh1_hav;sarm_veh2_hav;sarm_veh1_laf;sarm_veh2_laf; ...
    sarm_dcz1_hav;sarm_dcz2_hav;sarm_dcz1_laf;sarm_dcz2_laf];

sessd99drug_nocno = [sessd99veh;sessd99veh;sessd99veh;sessd99veh; ...
    sessd99dcz;sessd99dcz;sessd99dcz;sessd99dcz];
sessd99ses_nocno = [sessd99ses1;sessd99ses2;sessd99ses1;sessd99ses2; ...
    sessd99ses1;sessd99ses2;sessd99ses1;sessd99ses2];
sessd99subj_nocno = [sessd99hav;sessd99hav;sessd99laf;sessd99laf;sessd99hav; ...
    sessd99hav;sessd99laf;sessd99laf];
sessd99cnx_nocno = [sessd99connec;sessd99connec;sessd99connec;sessd99connec; ...
    sessd99connec;sessd99connec;sessd99connec;sessd99connec];

sesscharmdrug_nocno = [sesscharmveh;sesscharmveh;sesscharmveh;sesscharmveh; ...
    sesscharmdcz;sesscharmdcz;sesscharmdcz;sesscharmdcz];
sesscharmses_nocno = [sesscharmses1;sesscharmses2;sesscharmses1;sesscharmses2; ...
    sesscharmses1;sesscharmses2;sesscharmses1;sesscharmses2];
sesscharmsubj_nocno = [sesscharmhav;sesscharmhav;sesscharmlaf;sesscharmlaf;sesscharmhav; ...
    sesscharmhav;sesscharmlaf;sesscharmlaf];
sesscharmcnx_nocno = [sesscharmconnec;sesscharmconnec;sesscharmconnec;sesscharmconnec; ...
    sesscharmconnec;sesscharmconnec;sesscharmconnec;sesscharmconnec];

sesssarmdrug_nocno = [sesssarmveh;sesssarmveh;sesssarmveh;sesssarmveh; ...
    sesssarmdcz;sesssarmdcz;sesssarmdcz;sesssarmdcz];
sesssarmses_nocno = [sesssarmses1;sesssarmses2;sesssarmses1;sesssarmses2; ...
    sesssarmses1;sesssarmses2;sesssarmses1;sesssarmses2];
sesssarmsubj_nocno = [sesssarmhav;sesssarmhav;sesssarmlaf;sesssarmlaf;sesssarmhav; ...
    sesssarmhav;sesssarmlaf;sesssarmlaf];
sesssarmcnx_nocno = [sesssarmconnec;sesssarmconnec;sesssarmconnec;sesssarmconnec; ...
    sesssarmconnec;sesssarmconnec;sesssarmconnec;sesssarmconnec];

nocno_nesting = [0 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 0]; %nest session in subject

[p_sessd99_nocno,tbl_sessd99_nocno,stats_sessd99_nocno] = ...
    anovan(sessalld99_nocno,{sessd99drug_nocno, sessd99ses_nocno, ...
    sessd99subj_nocno, sessd99cnx_nocno},'model', ...
    [1 0 0 0;0 0 1 0;1 0 1 0],'random',[4],'nested',nocno_nesting, ...
    'varnames',{'Drug','Session','Subject','Area'});


[p_sesscharm_nocno,tbl_sesscharm_nocno,stats_sesscharm_nocno] = ...
    anovan(sessallcharm_nocno,{sesscharmdrug_nocno,sesscharmses_nocno, ...
    sesscharmsubj_nocno,sesscharmcnx_nocno},'model', ...
    [1 0 0 0;0 0 1 0;1 0 1 0],'random',[4],'nested',nocno_nesting, ...
    'varnames',{'Drug','Session','Subject','Area'});


[p_sesssarm_nocno,tbl_sesssarm_nocno,stats_sesssarm_nocno] = ...
    anovan(sessallsarm_nocno,{sesssarmdrug_nocno,sesssarmses_nocno, ...
    sesssarmsubj_nocno,sesssarmcnx_nocno},'model',... 
    [1 0 0 0;0 0 1 0;1 0 1 0],'random',[4],'nested',nocno_nesting, ...
    'varnames',{'Drug','Session','Subject','Area'});



%% Compare VEH and CNO treatments

sessalld99_nodcz = [d99_veh1_hav;d99_veh2_hav;d99_veh1_laf;d99_veh2_laf; ...
    d99_cno1_hav;d99_cno1_laf];
sessallcharm_nodcz = [charm_veh1_hav;charm_veh2_hav;charm_veh1_laf;charm_veh2_laf; ...
    charm_cno1_hav;charm_cno1_laf];
sessallsarm_nodcz = [sarm_veh1_hav;sarm_veh2_hav;sarm_veh1_laf;sarm_veh2_laf; ...
    sarm_cno1_hav;sarm_cno1_laf];

sessd99drug_nodcz = [sessd99veh;sessd99veh;sessd99veh;sessd99veh; ...
    sessd99cno;sessd99cno];
sessd99subj_nodcz = [sessd99hav;sessd99hav;sessd99laf;sessd99laf; ...
    sessd99hav;sessd99laf];
sessd99cnx_nodcz = [sessd99connec;sessd99connec;sessd99connec;sessd99connec; ...
    sessd99connec;sessd99connec];

sesscharmdrug_nodcz = [sesscharmveh;sesscharmveh;sesscharmveh;sesscharmveh; ...
    sesscharmcno;sesscharmcno];
sesscharmsubj_nodcz = [sesscharmhav;sesscharmhav;sesscharmlaf;sesscharmlaf; ...
    sesscharmhav;sesscharmlaf];
sesscharmcnx_nodcz = [sesscharmconnec;sesscharmconnec;sesscharmconnec;sesscharmconnec; ...
    sesscharmconnec;sesscharmconnec];

sesssarmdrug_nodcz = [sesssarmveh;sesssarmveh;sesssarmveh;sesssarmveh; ...
    sesssarmcno;sesssarmcno];
sesssarmsubj_nodcz = [sesssarmhav;sesssarmhav;sesssarmlaf;sesssarmlaf; ...
    sesssarmhav;sesssarmlaf];
sesssarmcnx_nodcz = [sesssarmconnec;sesssarmconnec;sesssarmconnec;sesssarmconnec; ...
    sesssarmconnec;sesssarmconnec];

[p_sessd99_nodcz,tbl_sessd99_nodcz,stats_sessd99_nodcz] = ...
    anovan(sessalld99_nodcz,{sessd99drug_nodcz,sessd99subj_nodcz,...
    sessd99cnx_nodcz},'model',[1 0 0;0 1 0;1 1 0],'random',...
    [3],'varnames',{'Drug','Subject','Area'});

[p_sesscharm_nodcz,tbl_sesscharm_nodcz,stats_sesscharm_nodcz] = ...
    anovan(sessallcharm_nodcz,{sesscharmdrug_nodcz,sesscharmsubj_nodcz,...
    sesscharmcnx_nodcz},'model',[1 0 0;0 1 0;1 1 0],'random',[3],...
    'varnames',{'Drug','Subject','Area'});

[p_sesssarm_nodcz,tbl_sesssarm_nodcz,stats_sesssarm_nodcz] = ...
    anovan(sessallsarm_nodcz,{sesssarmdrug_nodcz,sesssarmsubj_nodcz,...
    sesssarmcnx_nodcz},'model',[1 0 0;0 1 0;1 1 0],'random',[3],'varnames',...
    {'Drug','Subject','Area'});


%% Figures

%DCZ data
figure(400)
d99sess_nocno1 = d99_sessgroup(1:2,:);
d99sess_nocno2 = d99_sessgroup(4:5,:);
d99sess_nocno = vertcat(d99sess_nocno1, d99sess_nocno2);
d99nocno = bar(d99_by_drug_mean(1:2,:),0.9,'LineWidth',1.5); hold on;
d99nocno.FaceColor = 'flat';
d99nocno.CData(1,:) = [0.215686274509804	0.494117647058824	0.721568627450980];
d99nocno.CData(2,:) = [0.301960784313725	0.686274509803922	0.290196078431373];
for i=1:2
    plot([i,i],[d99_by_drug_mean(i)+d99_by_drug_sem(i) d99_by_drug_mean(i)-d99_by_drug_sem(i)],'k-');
end
d1 = plot(d99sess_nocno(1:2,:),'d','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
d2 = plot(d99sess_nocno(3:4,:),'+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
f400ax = gca;
f400ax.Box = 'off';
f400ax.Color = 'none';
f400ax.LineWidth = 1.5;
xlim([0.5 2.5]); ylim([-0.05 0.04]);
xticks(1:2); xticklabels({'VEH','DCZ'}); yticks([-0.05 -0.04 -0.03 -0.02...
    -0.01 0 0.01 0.02 0.03 0.04]);
title('D99');


figure(500);
charmsess_nocno1 = charm_sessgroup(1:2,:);
charmsess_nocno2 = charm_sessgroup(4:5,:);
charmsess_nocno = vertcat(charmsess_nocno1, charmsess_nocno2);
sarmsess_nocno1 = sarm_sessgroup(1:2,:);
sarmsess_nocno2 = sarm_sessgroup(4:5,:);
sarmsess_nocno = vertcat(sarmsess_nocno1, sarmsess_nocno2);

subplot(2,1,1); hold on;
charmnocno = bar(charm_by_drug_mean(1:2,:),0.9,'LineWidth',1.5); hold on;
charmnocno.FaceColor = 'flat';
charmnocno.CData(1,:) = [0.215686274509804	0.494117647058824	0.721568627450980];
charmnocno.CData(2,:) = [0.301960784313725	0.686274509803922	0.290196078431373];
for i=1:2
    plot([i,i],[charm_by_drug_mean(i)+charm_by_drug_sem(i) charm_by_drug_mean(i)-charm_by_drug_sem(i)],'k-');
end
plot(charmsess_nocno(1:2,:),'d','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(charmsess_nocno(3:4,:),'+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
f500ax1 = gca;
f500ax1.Box = 'off';
f500ax1.Color = 'none';
f500ax1.LineWidth = 1.5;
xlim([0.5 2.5]); ylim([-0.04 0.05]);
xticks(1:2); xticklabels({'VEH','DCZ'}); yticks([-0.04 -0.03 -0.02...
    -0.01 0 0.01 0.02 0.03 0.04 0.05]);
title('CHARM');

subplot(2,1,2); hold on;
sarmnocno = bar(sarm_by_drug_mean(1:2,:),0.9,'LineWidth',1.5); hold on;
sarmnocno.FaceColor = 'flat';
sarmnocno.CData(1,:) = [0.215686274509804	0.494117647058824	0.721568627450980];
sarmnocno.CData(2,:) = [0.301960784313725	0.686274509803922	0.290196078431373];
for i=1:2
    plot([i,i],[sarm_by_drug_mean(i)+sarm_by_drug_sem(i) sarm_by_drug_mean(i)-sarm_by_drug_sem(i)],'k-');
end
plot(sarmsess_nocno(1:2,:),'d','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(sarmsess_nocno(3:4,:),'+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
f500ax2 = gca;
f500ax2.Box = 'off';
f500ax2.Color = 'none';
f500ax2.LineWidth = 1.5;
xlim([0.5 2.5]); ylim([-0.02 0.07]);
xticks(1:2); xticklabels({'VEH','DCZ'}); yticks([-0.02...
    -0.01 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07])
title('SARM');




% CNO data
figure(401)
hav_d99sess_veh1 = d99_sessgroup(1,1);
hav_d99sess_veh2 = d99_sessgroup(1,2);
hav_d99sess_cno = [NaN d99_sessgroup(3,1)];
laf_d99sess_veh1 = d99_sessgroup(4,1);
laf_d99sess_veh2 = d99_sessgroup(4,2);
laf_d99sess_cno = [NaN d99_sessgroup(6,1)];
d99_drug_mean_nodcz = [d99_by_drug_mean(1,1);d99_by_drug_mean(3,1)];
d99_drug_sem_nodcz = [d99_by_drug_sem(1,1);d99_by_drug_sem(1,3)];
d99nodcz = bar(d99_drug_mean_nodcz(1:2,:),0.9,'LineWidth',1.5); hold on;
d99nodcz.FaceColor = 'flat';
d99nodcz.CData(1,:) = [0.215686274509804	0.494117647058824	0.721568627450980];
d99nodcz.CData(2,:) = [0.894117647058824	0.101960784313725	0.109803921568627];
for i=1:2
    plot([i,i],[d99_drug_mean_nodcz(i)+d99_drug_sem_nodcz(i) ...
        d99_drug_mean_nodcz(i)-d99_drug_sem_nodcz(i)],'k-');
end
plot(hav_d99sess_veh1(1,1),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(hav_d99sess_veh2(1,1),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(hav_d99sess_cno(1,:),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_d99sess_veh1(1,1),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_d99sess_veh2(1,1),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_d99sess_cno(1,:),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
f401ax = gca;
f401ax.Box = 'off';
f401ax.Color = 'none';
f401ax.LineWidth = 1.5;
xlim([0.5 2.5]);
ylim([-0.05 0.04]); xticks(1:2); xticklabels({'VEH','CNO'}); title('D99');


figure(501);
hav_charmsess_veh1 = charm_sessgroup(1,1);
hav_charmsess_veh2 = charm_sessgroup(1,2);
hav_charmsess_cno = [NaN charm_sessgroup(3,1)];
laf_charmsess_veh1 = charm_sessgroup(4,1);
laf_charmsess_veh2 = charm_sessgroup(4,2);
laf_charmsess_cno = [NaN charm_sessgroup(6,1)];
charm_drug_mean_nodcz = [charm_by_drug_mean(1,1);charm_by_drug_mean(3,1)];
charm_drug_sem_nodcz = [charm_by_drug_sem(1,1);charm_by_drug_sem(1,3)];
hav_sarmsess_veh1 = sarm_sessgroup(1,1);
hav_sarmsess_veh2 = sarm_sessgroup(1,2);
hav_sarmsess_cno = [NaN sarm_sessgroup(3,1)];
laf_sarmsess_veh1 = sarm_sessgroup(4,1);
laf_sarmsess_veh2 = sarm_sessgroup(4,2);
laf_sarmsess_cno = [NaN sarm_sessgroup(6,1)];
sarm_drug_mean_nodcz = [sarm_by_drug_mean(1,1);sarm_by_drug_mean(3,1)];
sarm_drug_sem_nodcz = [sarm_by_drug_sem(1,1);sarm_by_drug_sem(1,3)];

subplot(2,1,1); hold on;
charmnodcz = bar(charm_drug_mean_nodcz(1:2,:),0.9,'LineWidth',1.5); hold on;
charmnodcz.FaceColor = 'flat';
charmnodcz.CData(1,:) = [0.215686274509804	0.494117647058824	0.721568627450980];
charmnodcz.CData(2,:) = [0.894117647058824	0.101960784313725	0.109803921568627];
for i=1:2
    plot([i,i],[charm_drug_mean_nodcz(i)+charm_drug_sem_nodcz(i) ...
        charm_drug_mean_nodcz(i)-charm_drug_sem_nodcz(i)],'k-');
end
plot(hav_charmsess_veh1(1,1),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(hav_charmsess_veh2(1,1),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(hav_charmsess_cno(1,:),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_charmsess_veh1(1,1),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_charmsess_veh2(1,1),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_charmsess_cno(1,:),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
f501ax1 = gca;
f501ax1.Box = 'off';
f501ax1.Color = 'none';
f501ax1.LineWidth = 1.5;
xlim([0.5 2.5]);
ylim([-0.05 0.04]); xticks(1:2); xticklabels({'VEH','CNO'}); title('CHARM');

subplot(2,1,2); hold on;
sarmnodcz = bar(sarm_drug_mean_nodcz(1:2,:),0.9,'LineWidth',1.5); hold on;
sarmnodcz.FaceColor = 'flat';
sarmnodcz.CData(1,:) = [0.215686274509804	0.494117647058824	0.721568627450980];
sarmnodcz.CData(2,:) = [0.894117647058824	0.101960784313725	0.109803921568627];
for i=1:2
    plot([i,i],[sarm_drug_mean_nodcz(i)+sarm_drug_sem_nodcz(i) ...
        sarm_drug_mean_nodcz(i)-sarm_drug_sem_nodcz(i)],'k-');
end
plot(hav_sarmsess_veh1(1,1),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(hav_sarmsess_veh2(1,1),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(hav_sarmsess_cno(1,:),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_sarmsess_veh1(1,1),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_sarmsess_veh2(1,1),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(laf_sarmsess_cno(1,:),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
f501ax2 = gca;
f501ax2.Box = 'off';
f501ax2.Color = 'none';
f501ax2.LineWidth = 1.5;
xlim([0.5 2.5]);
ylim([-0.05 0.04]); xticks(1:2); xticklabels({'VEH','CNO'}); title('SARM');
