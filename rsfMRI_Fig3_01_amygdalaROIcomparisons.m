%Author: Catherine Elorette, Icahn School of Medicine at Mount Sinai
%This script requires the SPM software package
%This script uses amygdala ROI correlation files that have previously been generated in AFNI using standard commands (i.e. 3dTcorr1D)

clear all; close all;

%Get list of files
ROIdir = 'C:/Users/Catherine/Documents/Data/rsMRI/amyroicorrs/';
list = dir(ROIdir);
list = {list(:).name};

%Get the brain mask and convert to matrix format
maskname = list(contains(list,'mask'));
maskFile = char(maskname);
mask_pointer = spm_vol(fullfile(ROIdir,maskFile));
mask_vol = spm_read_vols(mask_pointer);
DataOrient = size(mask_vol);

mask_lin = mask_vol(:);
mask_loc = find(mask_lin==1);
mask_numvox = length(mask_loc);

%Prepare for ANOVAs later on
onevector = [];
onevector(1,1:mask_numvox) = 1;
twovector = [];
twovector(1,1:mask_numvox) = 2;

%Prepare file list
%Post DREADD surgery, already subtracted Post-Pre injection ('difference')
PDdiffiles = list(contains(list,'Dif.nii'));

%% Comparing the already-subtracted data

%Use the brain mask to isolate only voxels in the brain from the data files
for cc = 1:length(PDdiffiles)
    
    name = string(PDdiffiles(1,cc));
    funcFile = char(name);
    
    func_pointer = spm_vol(fullfile(ROIdir,funcFile));
    func_Data = spm_read_vols(func_pointer);
    
    func_Lin = reshape(func_Data,prod(DataOrient),size(func_Data,4));
    func_masked = func_Lin(mask_loc,:);
    
    difscans{cc} = func_masked;
end
% Order of difscans files to make analysis easier:
% CNO bilateral: 1, 16
% CNO Left: 4, 19
% CNO Right: 9, 24
% DCZ bilateral: 2, 3, 17, 18
% DCZ Left: 5, 6, 20, 21
% DCZ Right: 10, 11, 25, 26
% VEH bilateral: 14, 15, 29, 30
% VEH Left: 7, 8, 22, 23
% VEH Right: 12, 13, 27, 28


% Let's arrange our ANOVA by factors of: Drug, Session, Animal. So
% we'll start with Havoc DCZ 1, then Havoc DCZ 2, then Havoc VEH 1, etc.
% The drug order will therefore be DCZ then VEH, [1 1 2 2 1 1 2 2 ]. The session
% order will be [1 2 1 2 1 2 1 2]. The animal order will be Havoc = 1, 
% Lafeyette = 2, [1 1 1 1 2 2 2 2]. 

drugorder = [onevector onevector twovector twovector onevector onevector ...
    twovector twovector];
sessionorder = [onevector twovector onevector twovector onevector ...
    twovector onevector twovector];
animalorder = [onevector onevector onevector onevector ...
    twovector twovector twovector twovector];

difsession_nesting = [0 0 0; 0 0 1; 0 0 0]; %Nest session under animal

%Take the files where the bilateral amygdala ROI is used
difvector = [difscans{2}' difscans{3}' difscans{14}' difscans{15}' ...
    difscans{17}' difscans{18}' difscans{29}' difscans{30}'];


[dif_anovap, dif_anovatbl, dif_anovastats] = ...
    anovan(difvector,{drugorder sessionorder animalorder},'model', ...
    [1 0 0;0 0 1;1 0 1], 'nested', difsession_nesting, ...
    'varnames',{'Drug','Session','Animal'});

%Get the mean and sem values for each drug treatment so we can plot them
dczall = [difscans{2} difscans{3} difscans{17} difscans{18}];
dczallvector = [difscans{2}; difscans{3}; difscans{17}; difscans{18}];
difmean_dcz = mean(dczall,'all');
difmean_each_dcz = mean(dczall,1);
difsem_dcz = std(dczallvector,'omitnan')/(length(dczallvector).^0.5);

vehall = [difscans{14} difscans{15} difscans{29} difscans{30}];
vehallvector = [difscans{14}; difscans{15}; difscans{29}; difscans{30}];
difmean_veh = mean(vehall,'all');
difmean_each_veh = mean(vehall,1);
difsem_veh = std(vehallvector,'omitnan')/(length(vehallvector).^0.5);


dczvehvectors = [vehallvector'; dczallvector'];
difmeansdczveh = [difmean_veh; difmean_dcz];
difsemsdczveh = [difsem_veh; difsem_dcz];
havmeans = [difmean_each_veh(1) difmean_each_veh(2);difmean_each_dcz(1) difmean_each_dcz(2)];
lafmeans = [difmean_each_veh(3) difmean_each_veh(4);difmean_each_dcz(3) difmean_each_dcz(4)];

%Figure 3E
figure(100)
avvals2 = bar(difmeansdczveh(1:2,:),0.9,'LineWidth',1.5); hold on;
avvals2.FaceColor = 'flat';
avvals2.CData(1,:) = [0.215686274509804	0.494117647058824	0.721568627450980];
avvals2.CData(2,:) = [0.301960784313725	0.686274509803922	0.290196078431373];
for i=1:2
    plot([i,i],[difmeansdczveh(i)+difsemsdczveh(i) difmeansdczveh(i)-difsemsdczveh(i)],'k-','LineWidth',1.5);
end
plot(havmeans(1:2,:),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(lafmeans(1:2,:),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
f100ax = gca;
f100ax.Box = 'off';
f100ax.Color = 'none';
f100ax.LineWidth = 1.5;
xlim([0.5 2.5]);
ylim([-0.01 0.03]); yticks([-0.01 0 0.01 0.02 0.03]);
xticks(1:2); xticklabels({'VEH','DCZ'});


%% Compare CNO and VEH in the same way

% Let's arrange our second ANOVA by factors of: Drug, Session, Animal. So
% we'll start with Havoc CNO, then Havoc VEH 1, etc.
% The drug order will therefore be CNO then VEH, [1 2 2 1 2 2 ]. The session
% order will be [1 1 2 1 1 2]. The animal order will be Havoc = 1, 
% Lafeyette = 2, [1 1 1 2 2 2]. 

cnodrugorder = [onevector twovector twovector onevector ...
    twovector twovector];
cnosessionorder = [onevector onevector twovector onevector ...
    onevector twovector];
cnoanimalorder = [onevector onevector onevector ...
    twovector twovector twovector];

cnodifvector = [difscans{1}' difscans{14}' difscans{15}' ...
    difscans{16}' difscans{29}' difscans{30}'];

[cnodif_anovap, cnodif_anovatbl, cnodif_anovastats] = ...
    anovan(cnodifvector,{cnodrugorder cnosessionorder cnoanimalorder},'model', ...
    [1 0 0;0 0 1;1 0 1], 'nested', difsession_nesting, ...
    'varnames',{'Drug','Session','Animal'});

cnoall = [difscans{1} difscans{16}];
cnoallvector = [difscans{1}; difscans{16}];
difmean_cno = mean(cnoall,'all');
difmean_each_cno = mean(cnoall,1);
difsem_cno = std(cnoallvector,'omitnan')/(length(cnoallvector).^0.5);


difmeanscnoveh = [difmean_veh; difmean_cno];
difsemscnoveh = [difsem_veh; difsem_cno];
havmeanscno = [difmean_each_veh(1) difmean_each_veh(2);difmean_each_cno(1) NaN];
lafmeanscno = [difmean_each_veh(3) difmean_each_veh(4);difmean_each_cno(2) NaN];

%Figure 6B
figure(200)
cnoavvals2 = bar(difmeanscnoveh(1:2,:),0.9,'LineWidth',1.5); hold on;
cnoavvals2.FaceColor = 'flat';
cnoavvals2.CData(1,:) = [0.215686274509804	0.494117647058824	0.721568627450980];
cnoavvals2.CData(2,:) = [0.894117647058824	0.101960784313725	0.109803921568627];
for i=1:2
    plot([i,i],[difmeanscnoveh(i)+difsemscnoveh(i) difmeanscnoveh(i)-difsemscnoveh(i)],'k-','LineWidth',1.5);
end
plot(havmeanscno(1:2,:),'kd','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
plot(lafmeanscno(1:2,:),'k+','LineWidth',1.5,'MarkerSize',8,'Color',[0.3020 0.3020 0.3020]);
f100ax = gca;
f100ax.Box = 'off';
f100ax.Color = 'none';
f100ax.LineWidth = 1.5;
xlim([0.5 2.5]);
xticks(1:2); xticklabels({'VEH','CNO'});
ylim([-0.006 0.006]);
yticks([-0.006:0.002:0.006]);


