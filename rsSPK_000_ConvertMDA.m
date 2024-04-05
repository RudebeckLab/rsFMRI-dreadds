%% Converting plexon files to Matlab and MDA, needed for Spike sorting and further processing
%-
%- not needed to reproduce results, files are already processed
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023

clear

session = input('Session to convert into .mat files : ','s');
if strcmp(session(1),'G')
    cd(['I:\GENGER\' session]);

else
    cd(['I:\Havoc_REC\' session]);
end
%% Extracting Keyboard Inputs == Injection markers
if exist([session '_spl_KBD1.pl2'])==2
    disp('------------------------------------');
    disp('Extracting Events');
    StartingFileName = [session '_spl_KBD1.pl2'];
    readall
    EVT = [tsevs{36} ];
    save([OpenedFileName(1:end-4) '.mat'],'EVT','-v7.3')
    clear tsevs
end

if exist([session '_spl_KBD2.pl2'])==2
    disp('------------------------------------');
    disp('Extracting Events');
    StartingFileName = [session '_spl_KBD2.pl2'];
    readall
    EVT = [tsevs{37} ];
    save([OpenedFileName(1:end-4) '.mat'],'EVT','-v7.3')
    clear tsevs
end

if exist([session '_spl_KBD3.pl2'])==2
    disp('------------------------------------');
    disp('Extracting Events');
    StartingFileName = [session '_spl_KBD3.pl2'];
    readall
    EVT = [tsevs{38} ];
    save([OpenedFileName(1:end-4) '.mat'],'EVT','-v7.3')
    clear tsevs
end


%% Extracting Events from the Strobed.pl2 file
if exist([session '_spl_Strobed.mat'])~=2
    disp('------------------------------------');
    disp('Extracting Events');
    StartingFileName = [session '_spl_Strobed.pl2'];
    readall
    EVT = [tsevs{33} svStrobed];
    EVT(EVT(:,2)==0,:)=[];
    save([OpenedFileName(1:end-4) '.mat'],'EVT','-v7.3')

    clearvars -except elecs session
end

%% Extract EOG info from .pl2 files and save in .mat

list = dir([session '_spl_AI*.pl2']);

if ~isempty(list)
    disp('------------------------------------');
    disp('Extracting EOG');

    clear eog
    for i = 1 : 4 %length(list) %- only look for the first 5 analog channels (X, Y, pupil H, pupil W) + heart rate
        
        StartingFileName = list(i).name;
        disp('------------------------------------');
        disp(['Processing AI ' num2str(i) '/' num2str(5)]);
        readall
        
        clearvars -except allad adfreq OpenedFileName list f ich elecs i session m Row All_elecs all_data data eog
        
        keep = [];
        for j = 1 : length(allad)
            if ~isempty(allad{j})
                keep = [keep, j];
            end
        end
        
        eog(:,i) = allad{keep};
        clear allad
    end
    for i = 1 : length(list) %- delete AI .pl2 files
        delete(list(i).name)
    end
    eog = single(eog);
    Fs_eog = adfreq;
    save([session '_EOG.mat'],'eog','Fs_eog','-v7.3')
end


%% Converting .pl2 files in .MDA files (for MOUNTAINSORT in Ubuntu)
list = dir([session '_spl_WB*.pl2']);

for i = 1 : length(list)
    
    StartingFileName = list(i).name;
    disp('------------------------------------');
    disp(['Processing elec ' num2str(i) '/' num2str(length(list))]);
    readall
    
    clearvars -except allad adfreq OpenedFileName list f ich elecs i session m Row All_elecs all_data
    
    keep = [];
    for j = 1 : length(allad)
        if ~isempty(allad{j})
            keep = [keep, j];
        end
    end
    
    data = allad{keep};
    data = single(data);
    clear allad
    writemda16i(data',[OpenedFileName(1:end-4) '.mda']);
    delete(list(i).name)
    clear data
    
end
