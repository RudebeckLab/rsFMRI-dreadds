function [Data] = NCConnectomesGeneration(Atlas,level,subjname,diff_scaling)
% Comparing Connectomes 

if nargin < 2 
    disp('must enter Atlas and level variables')
    return
end


%level begins at 0 instead of 1 in the indexing for some reason, so level
%of '002' is actually level 3
if isempty(level)
    level = '002';
end

if length(level) ~= 3
    disp('Level needs to be string of 3 numbers')
    return
end

level = ['_' level '.'];


if isempty(Atlas)
    Atlas = 'D99';
end


%%
topdir = 'C:/Users/Catherine/Documents/Data/rsMRI/Connectomes/';
cd(topdir)
if string(Atlas) == 'CHARM'
    currdir = 'C:/Users/Catherine/Documents/Data/rsMRI/Connectomes/CHARM/';
elseif string(Atlas) == 'SARM'
    currdir = 'C:/Users/Catherine/Documents/Data/rsMRI/Connectomes/SARM/';
elseif string(Atlas) == 'D99'
    currdir = 'C:/Users/Catherine/Documents/Data/rsMRI/Connectomes/D99/';
end
cd(currdir)

%%levels = {'000';'001';'002';'003';'004';'005';};
conditions = {'PostDCZ1','PreDCZ1','PostDCZ2','PreDCZ2','PostVEH1','PreVEH1','PostVEH2','PreVEH2','PostCNO','PreCNO'};

if nargin < 3 || isempty(subjname) || strcmp(subjname,'all') || strcmp(subjname,'All')
    subj_list = dir(currdir);
    subj_list = {subj_list(3:end).name};
    subj_list = subj_list(~contains(subj_list,'DS'));
    subj_list = subj_list(~contains(subj_list,'Strip'));
    subj_list = subj_list(~contains(subj_list,'Results'));
else
    subj_list = {subjname};
end

filelist = dir([currdir subj_list{1} '/PreVEH1_' Atlas level 'netcc']);

tmp_file = [currdir '/' subj_list{1} '/' filelist(1).name];

tmp1 = fopen(tmp_file);
for ii = 1:5
    line = fgetl(tmp1);
end
C = textscan(line,'%d');
region_nums = double(C{1});
fclose(tmp1);
clear tmp1 C line


mn_r_vals = zeros(length(region_nums),length(region_nums),length(conditions),length(subj_list))*NaN;
mn_z_vals = mn_r_vals;
% (x,y,z,w) = regions, regions, conditions, subjects
clear tmp_file tmp
%% import data

for iii = 1:length(conditions)
    cur_cond = conditions{iii};

    for jjj = 1:length(subj_list)
        cur_subj = subj_list{jjj};
        cur_dir = [currdir cur_subj];
        cd(cur_dir)

        file_list = dir(cur_dir);
        file_list = {file_list(3:end).name}';
        net_files = contains(file_list,'.netcc');
        file_list = file_list(net_files);
        
        lvl_files = contains(file_list,level);
        lvl_files = file_list(lvl_files);

        cnd_files = contains(lvl_files,cur_cond);
        cnd_files = lvl_files(cnd_files);
        if ~isempty(cnd_files)
            for ijk = 1:length(cnd_files)
                tmp_file = cnd_files{ijk};
                disp(tmp_file)
                tmp = readmatrix(tmp_file,'FileType','text','NumHeaderLines',6);
                %start with the first file and use it as a template
                
                breaks = find(isnan(tmp(:,1)));
                %look for where there's a gap between the r values and z
                %scores

                mn_r_vals(:,:,iii,jjj) = tmp(1:breaks(1)-1,:);
                mn_z_vals(:,:,iii,jjj) = tmp(breaks(1)+1:end,:);
                %(region, region, conditions, subjects)
                %fill the r value matrix with everything above the gap, and
                %the z value matrix with everything below it
                %each loop fills up the matrix with values for each
                %condition and level
                clear tmp_file tmp breaks 
            end
            
        else
            disp([subj_list{jjj} ' ' cur_cond ' does not exist']) 
        end 
        
        cd(currdir)
        clear lvl_files cnd_files file_list cur_subj 
    end
    clear cur_cond 
end

clear iii ijk jjj cur_dir 


%% compare across DCZs
orig_scaling = [-1.25 1.25];

if nargin<4 || isempty(diff_scaling)
    diff_scaling = [-.2 .2];
end



for ijk = 1:length(subj_list)
    cur_mn_z = squeeze(mn_z_vals(:,:,:,ijk));
    diff_dcz1 = squeeze(cur_mn_z(:,:,1)) - squeeze(cur_mn_z(:,:,2));
    mn_dcz1(:,:,ijk) = diff_dcz1;
    diff_dcz2 = squeeze(cur_mn_z(:,:,3)) - squeeze(cur_mn_z(:,:,4));
    mn_dcz2(:,:,ijk) = diff_dcz2;
    diff_veh1 = squeeze(cur_mn_z(:,:,5)) - squeeze(cur_mn_z(:,:,6));
    mn_veh1(:,:,ijk) = diff_veh1;
    diff_veh2 = squeeze(cur_mn_z(:,:,7)) - squeeze(cur_mn_z(:,:,8));
    mn_veh2(:,:,ijk) = diff_veh2;
    diff_cno1 = squeeze(cur_mn_z(:,:,9)) - squeeze(cur_mn_z(:,:,10));
    mn_cno1(:,:,ijk) = diff_cno1;
end
clear ijk 
%Subtract Pre injection z and r values from Post injection z and r values,
%then store the results in a 3-D matrix where (connection, connection,
%subject) is the format. 


%% anova
drug_var = [];
subj_var = [];
area_var = [];
sess_var = [];
%make sure to use only half the matrix 
data = [];
for ijk = 1:length(subj_list)
    tmp_matrix = tril(ones(size(mn_veh1,1),size(mn_veh1,2)),-1);
    
    % get subject data
    cur_vehl = squeeze(mn_veh1(:,:,ijk));
    cur_veh2 = squeeze(mn_veh2(:,:,ijk));
    cur_dcz1 = squeeze(mn_dcz1(:,:,ijk));
    cur_dcz2 = squeeze(mn_dcz2(:,:,ijk));
    cur_cno1 = squeeze(mn_cno1(:,:,ijk));
    % use only one side of the matrix
    cur_vehl = tril(cur_vehl);
    cur_dcz1 = tril(cur_dcz1);
    cur_veh2 = tril(cur_veh2);
    cur_dcz2 = tril(cur_dcz2);
    cur_cno1 = tril(cur_cno1);
    % nan out the extra values
    cur_vehl(tmp_matrix==0)= NaN;
    cur_veh2(tmp_matrix==0)= NaN;
    cur_dcz1(tmp_matrix==0)= NaN;
    cur_dcz2(tmp_matrix==0)= NaN;
    cur_cno1(tmp_matrix==0)= NaN;
    % turn it into a linear array for concantenation 
    cur_vehl = cur_vehl(:);
    cur_veh2 = cur_veh2(:);
    cur_dcz1 = cur_dcz1(:);
    cur_dcz2 = cur_dcz2(:);
    cur_cno1 = cur_cno1(:);
    % get rid of the nans 
    cur_vehl = cur_vehl(~isnan(cur_vehl));
    cur_veh2 = cur_veh2(~isnan(cur_veh2));
    cur_dcz1 = cur_dcz1(~isnan(cur_dcz1));
    cur_dcz2 = cur_dcz2(~isnan(cur_dcz2));
    cur_cno1 = cur_cno1(~isnan(cur_cno1));
    disp(ijk) 
    disp([length(cur_vehl) length(cur_veh2) length(cur_dcz1) length(cur_dcz2) length(cur_cno1)])
    
    data = [data; cur_vehl; cur_veh2; cur_dcz1; cur_dcz2; cur_cno1;];
    drug_var = [drug_var; ones(size(cur_vehl)); ones(size(cur_veh2)); ones(size(cur_dcz1))*2; ones(size(cur_dcz2))*2; ones(size(cur_cno1))*3];
    subj_var = [subj_var; ones(length(cur_vehl)+length(cur_veh2)+length(cur_dcz1)+length(cur_dcz2)+length(cur_cno1),1)*ijk];
    tmp_area = [1:length(cur_vehl) 1:length(cur_veh2) 1:length(cur_dcz1) 1:length(cur_dcz2) 1:length(cur_cno1)]';
    area_var = [area_var; tmp_area];
    sess_var = [sess_var; ones(size(cur_vehl))*1; ones(size(cur_veh2))*2; ones(size(cur_dcz1))*1; ones(size(cur_dcz2))*2; ones(size(cur_cno1))*1];
    clear cur* tmp_area 
end

cd('C:/Users/Catherine/Documents/Data/rsMRI/')

Data = {'Zvalues', 'Subject', 'Drug', 'Session', 'Connection', 'Region_Nums';
    data, subj_var, drug_var, sess_var, area_var, region_nums};
if string(Atlas) == 'CHARM'
    save('CHARMData.mat','Data');
elseif string(Atlas) == 'SARM'
    save('SARMData.mat','Data');
elseif string(Atlas) == 'D99'
    save('D99Data.mat','Data');
end
