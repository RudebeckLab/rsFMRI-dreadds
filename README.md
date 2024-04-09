
# rsFMRI-dreadds

**Resting State fMRI analyses:**

These scripts are meant to be used with preprocessed data derived from standard AFNI pipelines. The dates and order of experiments are as follows: 2020/03/16 Animal L DCZ 1; 2020/05/27 Animal H vehicle 1; 2020/06/01 Animal L vehicle 1; 2020/06/08 Animal H DCZ 1; 2020/06/22 Animal L DCZ 2; 2020/06/24 Animal H vehicle 2; 2020/07/06 Animal L vehicle 2; 2020/07/08 Animal H DCZ 2; 2021/10/25 Animal L CNO; 2022/01/13 Animal H CNO.

For atlas connectome analysis, related to Figure 4A-C and Figure 6C: files are organized in folders according to the atlas used, then in folders by subject. File names should be organized as 'Pre/Post' 'Treatment' '1/2' _ 'Atlas' _ 'Level' .netcc. 

Requirement: 
  -  Standard atlases in NMT space, available from PRIME-DRE: https://prime-re.github.io/templates_and_atlases.html#macaque_atlases

For amygdala ROI correlation anaylses, related to Figure 3E and Figure 6B: this script uses files for amygdala ROI correlation maps generated in AFNI, then corrected for within-session comparisons (subtract the pre-injection correlation map from the post-injection correlation map). Files should be organized as 'Subject' _ 'Treatment' _ '1/2' _ 'Dif' .nii

Requirement:
  -  SPM software package: https://www.fil.ion.ucl.ac.uk/spm/

For voxel to voxel correlation analyses, related to Figure 5E: this script uses ROIs defined from standard atlases in NMT space. 

Requirement:
  -  SPM software package: https://www.fil.ion.ucl.ac.uk/spm/
  -  Standard atlases in NMT space, available from PRIME-DRE: https://prime-re.github.io/templates_and_atlases.html#macaque_atlases


**Resting State SPK and LFP analyses:**

The dates and order of experiments are as follows: 2022/04/14 Animal H DCZ 1; 2022/04/21 Animal H vehicle 1; 2022/05/12 Animal H DCZ 2; 2022/05/24 Animal H vehicle 2; 2022/06/03 Animal H CNO; 2022/08/18 Animal L vehicle 1; 2022/08/25 Animal L DCZ 1; 2022/09/01 Animal L DCZ 2; 2022/09/08 Animal L vehicle 2; 2022/09/22 Animal L CNO.

Run scripts in order 001 --> 002 ...
Scripts "..._000_ConvertMDA.m" and "..._000_ExtractSpikes.m" are not required if starting from sorted neurons and filtered LFPs. 

Requirement : 

- Fieldtrip : https://github.com/fieldtrip/fieldtrip

- Spike and LFP files!

For running sorting and initial conversion steps :

- Plexon Offline Files SDK (for reading pl2 files) : https://plexon.com/software-downloads/#software-downloads-SDKs

- MountainSort input/output matlab functions (readmda..) : https://github.com/flatironinstitute/mountainlab-js/tree/master/utilities/matlab

Additional tools are in the "utils" folder. Some functions there were not created by myself and might be under some Copyright. See individual functions for that information.

