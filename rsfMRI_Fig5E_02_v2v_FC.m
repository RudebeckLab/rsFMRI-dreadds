function [corr_ROI1_ROI2] = rsfMRI_Fig5E_02_v2v_FC(funcFile,ROI1file,ROI2file)

% ROI-to-ROI voxel analysis
ROIdir = '/Users/Catherine/Documents/Data/rsMRI/ROIs/HAVLAF/';

%% Load ROIs

ROI1_pointer = spm_vol(fullfile(ROIdir,ROI1file));
ROI1_vol = spm_read_vols(ROI1_pointer);

ROI2_pointer = spm_vol(fullfile(ROIdir,ROI2file));
ROI2_vol = spm_read_vols(ROI2_pointer);

DataOrient = size(ROI1_vol);


%% Load Rest Data

func_pointer = spm_vol(fullfile(ROIdir,funcFile));
func_Data = spm_read_vols(func_pointer);

func_Lin = reshape(func_Data,prod(DataOrient),size(func_Data,4)); 

%% Extract ROI data
% Amygdala
ROI1_lin = ROI1_vol(:);
ROI2_lin = ROI2_vol(:);

ROI1_loc = find(ROI1_lin==1); % gets index of the ROI voxels 
ROI2_loc = find(ROI2_lin==1); % gets index of the ROI voxels 

ROI1_numvox = length(ROI1_loc);
ROI2_numvox = length(ROI2_loc);

func_ROI1 = func_Lin(ROI1_loc,:);
func_ROI2 = func_Lin(ROI2_loc,:);


%% Calculate Voxel to Voxel Correlation 


all_func_data = [func_ROI1; func_ROI2];
all_func_flipped = all_func_data';
all_func_corr = corrcoef(all_func_flipped);

[all_func_corr,p_func_corr,RL_func_corr,RU_func_corr] = corrcoef(all_func_flipped);
rval_ROI1_ROI2 = all_func_corr(ROI1_numvox+1:(ROI1_numvox+ROI2_numvox),1:ROI1_numvox);
pval_ROI1_ROI2 = p_func_corr(ROI1_numvox+1:(ROI1_numvox+ROI2_numvox),1:ROI1_numvox);
RL_ROI1_ROI2 = RL_func_corr(ROI1_numvox+1:(ROI1_numvox+ROI2_numvox),1:ROI1_numvox);
RU_ROI1_ROI2 = RU_func_corr(ROI1_numvox+1:(ROI1_numvox+ROI2_numvox),1:ROI1_numvox);

pval05_ROI1_ROI2 = pval_ROI1_ROI2 < 0.05;
pvaltotal = pval_ROI1_ROI2(pval05_ROI1_ROI2);


all_func_zs = atanh(all_func_corr); 
%makes the rows the time series per voxel (columns). 
%corrcoef calculates correlation coefficients in which each row is an 
%observation and each column is a variable

%% Extract ROI to other 

% Extract correlations between ROIs Correlations 
corr_ROI1_ROI2 = all_func_zs(ROI1_numvox+1:(ROI1_numvox+ROI2_numvox),1:ROI1_numvox);



%% Calculate Mean Correlations

figure(210); clf;
imagesc(corr_ROI1_ROI2)
set(gca,'xtick',[],'ytick',[])
