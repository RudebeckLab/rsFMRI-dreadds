**testing**

# rsSPK-dreadds
Resting State SPK and LFP analyses

Run scripts in order 001 --> 002 ...
Scripts "..._000_ConvertMDA.m" and "..._000_ExtractSpikes.m" are not required if starting from sorted neurons and filtered LFPs. 

Requirement : 

- Fieldtrip : https://github.com/fieldtrip/fieldtrip

- Spike and LFP files!

For running sorting and initial conversion steps :

- Plexon Offline Files SDK (for reading pl2 files) : https://plexon.com/software-downloads/#software-downloads-SDKs

- MountainSort input/output matlab functions (readmda..) : https://github.com/flatironinstitute/mountainlab-js/tree/master/utilities/matlab

Additional tools are in the "utils" folder. Some functions there were not created by myself and might be under some Copyright. See individual functions for that information.

