%Script written by Brian Russ, revised by Catherine Elorette. Analysis was
%created by Vincent Costa and developed for our study with his advice and
%permission.

clear all; close all;

cd('/Users/Catherine/Documents/Data/rsMRI/ROIs/HAVLAF/');

Left_Amygdala_File = 'D99_LAMY.nii';
Right_Amygdala_File = 'D99_RAMY.nii';

Left_vlPFC_File = 'CHARM6_L12OL15.nii';
Right_vlPFC_File = 'CHARM6_R12OL15.nii';

%Change this to 'Left_Amygdala_File' and 'Left_vlPFC_File' to run the left
%side analysis
ROIs = {Right_Amygdala_File; Right_vlPFC_File};

%Load in the whole brain raw data files for each drug treatment before and
%after injection
HAVVEH1Pre_File = 'HAV20200527Pre.nii';
HAVVEH1Post_File = 'HAV20200527Post.nii';
HAVVEH2Pre_File = 'HAV20200624Pre.nii';
HAVVEH2Post_File = 'HAV20200624Post.nii';
HAVDCZ1Pre_File = 'HAV20200608Pre.nii';
HAVDCZ1Post_File = 'HAV20200608Post.nii';
HAVDCZ2Pre_File = 'HAV20200708Pre.nii';
HAVDCZ2Post_File = 'HAV20200708Post.nii';

LAFVEH1Pre_File = 'LAF20200601Pre.nii';
LAFVEH1Post_File = 'LAF20200601Post.nii';
LAFVEH2Pre_File = 'LAF20200706Pre.nii';
LAFVEH2Post_File = 'LAF20200706Post.nii';
LAFDCZ1Pre_File = 'LAF20200316Pre.nii';
LAFDCZ1Post_File = 'LAF20200316Post.nii';
LAFDCZ2Pre_File = 'LAF20200622Pre.nii';
LAFDCZ2Post_File = 'LAF20200622Post.nii';

%Organize the pre and post injection files for later subtraction
Pre_Errts_File = {HAVVEH1Pre_File; HAVVEH2Pre_File; HAVDCZ1Pre_File; HAVDCZ2Pre_File; LAFVEH1Pre_File; LAFVEH2Pre_File; LAFDCZ1Pre_File; LAFDCZ2Pre_File};
Post_Errts_File = {HAVVEH1Post_File; HAVVEH2Post_File; HAVDCZ1Post_File; HAVDCZ2Post_File; LAFVEH1Post_File; LAFVEH2Post_File; LAFDCZ1Post_File; LAFDCZ2Post_File};

Diff_corr_files = {'HAVVEH1'; 'HAVVEH2'; 'HAVDCZ1'; 'HAVDCZ2'; 'LAFVEH1'; 'LAFVEH2'; 'LAFDCZ1'; 'LAFDCZ2'};

%% Correlate each ROI with the other ROIs
%This section generates correlations between every voxel of ROI 1 and every
%voxel of ROI 2. Since the DREADDs were expressed in the amygdala, we keep
%amygdala as ROI 1. This allows us to look at changes in correlation from
%the perspective of the amygdala

Diff_corrs = {};
Pre_corrs_files = {};
Post_corrs_files = {};

strt_indx = 1;
for abc = 1:length(Post_Errts_File)
    for ijk = 1:length(ROIs)-1
        ROI1 = ROIs{ijk};
        for lmn = ijk+1:length(ROIs)
            ROI2 = ROIs{lmn};
            
            Corr_Order(:,strt_indx) = {ROI1; ROI2};
            Corr_Order(:,strt_indx+1) = {ROI2; ROI1};
            
            
            pre_corrs = rsfMRI_Fig5E_02_v2v_FC(Pre_Errts_File{abc},ROI1,ROI2);
            post_corrs = rsfMRI_Fig5E_02_v2v_FC(Post_Errts_File{abc},ROI1,ROI2);
            
            Pre_corrs_files{abc} = pre_corrs(:,:);
            Post_corrs_files{abc} = post_corrs(:,:);
            Diff_corrs{abc} = post_corrs(:,:) - pre_corrs(:,:);
            
            
            clear ROI2
            
            strt_indx = strt_indx+2;
            close all
        end
        clear ROI1 lmn
    end
    
end
clear strt_indx ijk 

%% Calculations


%You can check the order of the cells by looking at the Diff_corr_files
%array. If you want to analyze just the pre and post injection data, those
%get saved in the 'Pre_corrs_files' and 'Post_corrs_files' arrays
HAVVEH1 = Diff_corrs{1,1};
HAVVEH2 = Diff_corrs{1,2};
HAVDCZ1 = Diff_corrs{1,3};
HAVDCZ2 = Diff_corrs{1,4};
LAFVEH1 = Diff_corrs{1,5};
LAFVEH2 = Diff_corrs{1,6};
LAFDCZ1 = Diff_corrs{1,7};
LAFDCZ2 = Diff_corrs{1,8};

%This averages the changes due to drug injection over both sessions for
%each animal
HAVVEH_av = (Diff_corrs{1,1}+Diff_corrs{1,2})/2;
HAVDCZ_av = (Diff_corrs{1,3}+Diff_corrs{1,4})/2;
LAFVEH_av = (Diff_corrs{1,5}+Diff_corrs{1,6})/2;
LAFDCZ_av = (Diff_corrs{1,7}+Diff_corrs{1,8})/2;

%Now look at the effect of DCZ compared to vehicle
HAV_av = HAVDCZ_av - HAVVEH_av;
LAF_av = LAFDCZ_av - LAFVEH_av;

HAV_mean = mean(HAV_av,1);
LAF_mean = mean(LAF_av,1);

Both_mean = (HAV_mean+LAF_mean)/2;

%Separate the effect of DCZ vs vehicle by voxel so we see which voxels had
%an increase in correlation between our ROIs and which had a decrease
HAV_pos = HAV_av > 0.1;
HAV_neg = HAV_av < -0.1;
LAF_pos = LAF_av > 0.1;
LAF_neg = LAF_av < -0.1;

HAVpos_all = sum(HAV_pos);
HAVneg_all = sum(HAV_neg);
LAFpos_all = sum(LAF_pos);
LAFneg_all = sum(LAF_neg);
HAVpos_total = sum(HAVpos_all,2);
HAVneg_total = sum(HAVneg_all,2);
LAFpos_total = sum(LAFpos_all,2);
LAFneg_total = sum(LAFneg_all,2);

%Prepare labels
for x = 1:length(ROIs)
    Corr_List(x) = {strcat(Corr_Order(1,x),'_to_',Corr_Order(2,x))};
end

disp('Finished ROI to ROI correlations')

%% Visualize the correlation maps
%Optional

figure(901);
subplot(2,2,1);
imagesc(Diff_corrs{1,1});
title('HAV VEH1');
ax9011 = gca;
ax9011.YDir = 'normal'
ax9011.YLabel.String = 'ROI 2 voxels';
ax9011.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
subplot(2,2,2);
imagesc(Diff_corrs{1,2});
title('HAV VEH2');
ax9012 = gca;
ax9012.YDir = 'normal'
ax9012.YLabel.String = 'ROI 2 voxels';
ax9012.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
subplot(2,2,3);
imagesc(Diff_corrs{1,5});
title('LAF VEH1');
ax9013 = gca;
ax9013.YDir = 'normal'
ax9013.YLabel.String = 'ROI 2 voxels';
ax9013.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
subplot(2,2,4);
imagesc(Diff_corrs{1,6});
title('LAF VEH2');
ax9014 = gca;
ax9014.YDir = 'normal'
ax9014.YLabel.String = 'ROI 2 voxels';
ax9014.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);

figure(902);
subplot(2,2,1);
imagesc(Diff_corrs{1,3});
ax9021 = gca;
ax9021.YDir = 'normal'
ax9021.YLabel.String = 'ROI 2 voxels';
ax9021.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
title('HAV DCZ1');
subplot(2,2,2);
imagesc(Diff_corrs{1,4});
ax9022 = gca;
ax9022.YDir = 'normal'
ax9022.YLabel.String = 'ROI 2 voxels';
ax9022.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
title('HAV DCZ2');
subplot(2,2,3);
imagesc(Diff_corrs{1,7});
ax9023 = gca;
ax9023.YDir = 'normal'
ax9023.YLabel.String = 'ROI 2 voxels';
ax9023.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
title('LAF DCZ1');
subplot(2,2,4);
imagesc(Diff_corrs{1,8});
title('LAF DCZ2');
ax9024 = gca;
ax9024.YDir = 'normal'
ax9024.YLabel.String = 'ROI 2 voxels';
ax9024.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);

figure(903);
subplot(2,2,1);
imagesc(HAVVEH_av);
title('HAV VEH average');
ax9031 = gca;
ax9031.YDir = 'normal'
ax9031.YLabel.String = 'ROI 2 voxels';
ax9031.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
subplot(2,2,2);
imagesc(LAFVEH_av);
title('LAF VEH average');
ax9032 = gca;
ax9032.YDir = 'normal'
ax9032.YLabel.String = 'ROI 2 voxels';
ax9032.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
subplot(2,2,3);
imagesc(HAVDCZ_av);
title('HAV DCZ average');
ax9033 = gca;
ax9033.YDir = 'normal'
ax9033.YLabel.String = 'ROI 2 voxels';
ax9033.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);
subplot(2,2,4);
imagesc(LAFDCZ_av);
title('LAF DCZ average');
ax9034 = gca;
ax9034.YDir = 'normal'
ax9034.YLabel.String = 'ROI 2 voxels';
ax9034.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);

figure(904);
imagesc(HAV_av);
title('HAV DCZ-VEH average');
ax904 = gca;
ax904.YDir = 'normal'
ax904.YLabel.String = 'ROI 2 voxels';
ax904.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);

figure(905);
imagesc(LAF_av);
title('LAF DCZ-VEH average');
ax905 = gca;
ax905.YDir = 'normal'
ax905.YLabel.String = 'ROI 2 voxels';
ax905.XLabel.String = 'ROI 1 voxels';
colorbar('Ticks',[-0.2,-0.1,0,0.1,0.2],'Limits',[-0.2 0.2]);

%% Visualize results on the amygdala ROI
% Figure 5E
%Get a 3D amygdala mask and transform it into 2D space
ROIdir = '/Users/Catherine/Documents/Data/rsMRI/ROIs/HAVLAF/';
if string(ROIs(1)) == 'D99_RAMY.nii'
    AMY_pointer = spm_vol(fullfile(ROIdir,Right_Amygdala_File));
elseif string(ROIs(1)) == 'D99_LAMY.nii'
    AMY_pointer = spm_vol(fullfile(ROIdir,Left_Amygdala_File));
end
AMY_vol = spm_read_vols(AMY_pointer);
AMY_lin = AMY_vol(:);
AMY_loc = find(AMY_lin==1);

%Arrange the 2D correlation map of the amygdala back into 3D space
AMY_HAVmean = AMY_lin;
AMY_HAVmean(AMY_loc) = HAV_mean';
AMY_HAVmean3D = reshape(AMY_HAVmean,43,52,33);

AMY_LAFmean = AMY_lin;
AMY_LAFmean(AMY_loc) = LAF_mean';
AMY_LAFmean3D = reshape(AMY_LAFmean,43,52,33);

%There's not a great way to look at the 3D data, so let's look at the data
%by slicing it coronally
AMY_HAVmean3D_30slice = squeeze(AMY_HAVmean3D(:,30,:));
AMY_HAVmean3D_31slice = squeeze(AMY_HAVmean3D(:,31,:));
AMY_HAVmean3D_32slice = squeeze(AMY_HAVmean3D(:,32,:));
AMY_HAVmean3D_33slice = squeeze(AMY_HAVmean3D(:,33,:));
AMY_HAVmean3D_34slice = squeeze(AMY_HAVmean3D(:,34,:));

AMY_LAFmean3D_30slice = squeeze(AMY_LAFmean3D(:,30,:));
AMY_LAFmean3D_31slice = squeeze(AMY_LAFmean3D(:,31,:));
AMY_LAFmean3D_32slice = squeeze(AMY_LAFmean3D(:,32,:));
AMY_LAFmean3D_33slice = squeeze(AMY_LAFmean3D(:,33,:));
AMY_LAFmean3D_34slice = squeeze(AMY_LAFmean3D(:,34,:));

%Prepare just the amygdala mask 
AMY_vol_30slice = squeeze(AMY_vol(:,30,:));
AMY_vol_31slice = squeeze(AMY_vol(:,31,:));
AMY_vol_32slice = squeeze(AMY_vol(:,32,:));
AMY_vol_33slice = squeeze(AMY_vol(:,33,:));
AMY_vol_34slice = squeeze(AMY_vol(:,34,:));

%Visualize the data. Leftmost figures are the most posterior, rightmost the
%most anterior. The ROI is displayed rotated 90 degrees such that the top
%of the graph is the left side of the brain and the bottom of the graph is
%the right side.
figure(503);
tiledlayout(2,5)
nexttile
imagesc(AMY_HAVmean3D_30slice(:,:))
axis([0 20 12 32])
title('Havoc mean change slice 30')
nexttile
imagesc(AMY_HAVmean3D_31slice(:,:))
axis([0 20 12 32])
title('Havoc mean change slice 31')
nexttile
imagesc(AMY_HAVmean3D_32slice(:,:))
axis([0 20 12 32])
title('Havoc mean change slice 32')
nexttile
imagesc(AMY_HAVmean3D_33slice(:,:))
axis([0 20 12 32])
title('Havoc mean change slice 33')
nexttile
imagesc(AMY_HAVmean3D_34slice(:,:))
axis([0 20 12 32])
title('Havoc mean change slice 34')

nexttile
imagesc(AMY_LAFmean3D_30slice(:,:))
axis([0 20 12 32])
title('Lafeyette mean change slice 30')
nexttile
imagesc(AMY_LAFmean3D_31slice(:,:))
axis([0 20 12 32])
title('Lafeyette mean change slice 31')
nexttile
imagesc(AMY_LAFmean3D_32slice(:,:))
axis([0 20 12 32])
title('Lafeyette mean change slice 32')
nexttile
imagesc(AMY_LAFmean3D_33slice(:,:))
axis([0 20 12 32])
title('Lafeyette mean change slice 33')
nexttile
imagesc(AMY_LAFmean3D_34slice(:,:))
axis([0 20 12 32])
title('Lafeyette mean change slice 34')
colormap(parula(30));
clim([-0.05 0.05])
cb = colorbar('Ticks',[-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05],'Limits',[-0.05 0.05]);
cb.Layout.Tile = 'east';


figure(403);
tiledlayout(2,5)
nexttile
imagesc(AMY_vol_30slice(:,:))
axis([0 20 12 32])
title('Slice 30')
nexttile
imagesc(AMY_vol_31slice(:,:))
axis([0 20 12 32])
title('Slice 31')
nexttile
imagesc(AMY_vol_32slice(:,:))
axis([0 20 12 32])
title('Slice 32')
nexttile
imagesc(AMY_vol_33slice(:,:))
axis([0 20 12 32])
title('Slice 33')
nexttile
imagesc(AMY_vol_34slice(:,:))
axis([0 20 12 32])
title('Slice 34')

nexttile
imagesc(AMY_vol_30slice(:,:))
axis([0 20 12 32])
title('Slice 30')
nexttile
imagesc(AMY_vol_31slice(:,:))
axis([0 20 12 32])
title('Slice 31')
nexttile
imagesc(AMY_vol_32slice(:,:))
axis([0 20 12 32])
title('Slice 32')
nexttile
imagesc(AMY_vol_33slice(:,:))
axis([0 20 12 32])
title('Slice 33')
nexttile
imagesc(AMY_vol_34slice(:,:))
axis([0 20 12 32])
title('Slice 34')
cb = colorbar('Ticks',[-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05],'Limits',[-0.05 0.05]);
cb.Layout.Tile = 'east';
