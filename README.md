
# rsFMRI-dreadds

**Resting State fMRI analyses:**

These scripts are meant to be used with preprocessed data derived from standard AFNI pipelines. 

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

Run scripts in order 001 --> 002 ...
Scripts "..._000_ConvertMDA.m" and "..._000_ExtractSpikes.m" are not required if starting from sorted neurons and filtered LFPs. 

Requirement : 

- Fieldtrip : https://github.com/fieldtrip/fieldtrip

- Spike and LFP files!

For running sorting and initial conversion steps :

- Plexon Offline Files SDK (for reading pl2 files) : https://plexon.com/software-downloads/#software-downloads-SDKs

- MountainSort input/output matlab functions (readmda..) : https://github.com/flatironinstitute/mountainlab-js/tree/master/utilities/matlab

Additional tools are in the "utils" folder. Some functions there were not created by myself and might be under some Copyright. See individual functions for that information.

